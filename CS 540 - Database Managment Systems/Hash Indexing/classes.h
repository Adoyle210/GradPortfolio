#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <cmath>

using namespace std;

class Record {
public:
    int id, manager_id; // Employee ID and their manager's ID
    string bio, name; // Fixed length string to store employee name and biography

    Record(int id, string name, string bio, int manager_id)
        : id(id), name(name), bio(bio), manager_id(manager_id) {}

    // Constructor accepting a vector<string>
    Record(vector<string> &fields) {
        id = stoi(fields[0]);
        name = fields[1];
        bio = fields[2];
        manager_id = stoi(fields[3]);
    }

    // Function to get the size of the record
    int get_size() {
        return int(sizeof(int) * 4 + bio.length() + name.length());
    }

    // Function to serialize the record for writing to file
    string serialize() const {
        ostringstream oss;
        oss.write(reinterpret_cast<const char *>(&id), sizeof(id));
        oss.write(reinterpret_cast<const char *>(&manager_id), sizeof(manager_id));
        int name_len = name.size();
        int bio_len = bio.size();
        oss.write(reinterpret_cast<const char *>(&name_len), sizeof(name_len));
        oss.write(name.c_str(), name.size());
        oss.write(reinterpret_cast<const char *>(&bio_len), sizeof(bio_len));
        oss.write(bio.c_str(), bio.size());
        return oss.str();
    }

    void print() {
        cout << "\tID: " << id << "\n";
        cout << "\tNAME: " << name << "\n";
        cout << "\tBIO: " << bio << "\n";
        cout << "\tMANAGER_ID: " << manager_id << "\n";
    }
};

class Page {
public:
    vector<Record> records; // Data_Area containing the records
    vector<pair<int, int> > slot_directory; // Slot directory containing offset and size of each record
    int cur_size = 0; // Current size of the page
    int records_size = 0;
    int overflowPointerIndex; // Offset of overflow page, set to -1 by default

    // Constructor
    Page() : overflowPointerIndex(-1) {}

    void print(){
        cout << "# Records: " << records.size() << endl;
        for (auto& r: records) {
            r.print();
        }
    }

    // Function to insert a record into the page
    bool insert_record_into_page(Record r) {
        // Check if page size limit exceeded
        if (cur_size + r.get_size() + (sizeof(int) * 2) >= 4096) { 
            return false; 
        } else {
            records.push_back(r);
            slot_directory.push_back(make_pair(records_size, r.get_size()));
            records_size += r.get_size();
            cur_size += r.get_size() + (sizeof(int) * 2);
            return true;
        }
    }

    // Function to write the page to a binary output stream. You may use
    void write_into_data_file(ostream &out) const {
        char page_data[4096] = {0}; // Buffer to hold page data
        int offset = 0;
        
        // Write records into page_data buffer
        for (const auto &record: records) {
            string serialized = record.serialize();
            memcpy(page_data + offset, serialized.c_str(), serialized.size());
            offset += serialized.size();
        }

        //DONE:
        // Write the number of records contained in the page
        offset = 0;
        int num_records = slot_directory.size();
        memcpy(page_data + 4096 - sizeof(num_records), &num_records, sizeof(num_records));
        offset += sizeof(num_records);

        // Write overflowPointerIndex into page_data buffer
        memcpy(page_data + 4096 - sizeof(overflowPointerIndex) - offset, &overflowPointerIndex, sizeof(overflowPointerIndex));
        offset += sizeof(overflowPointerIndex);

        //  You should write the first entry of the slot_directory, which have the info about the first record at the bottom of the page, before overflowPointerIndex.
        for (const auto& slot : slot_directory) {
            memcpy(page_data + 4096 - offset - sizeof(slot.first), &slot.first, sizeof(slot.first));
            offset += sizeof(slot.first);

            memcpy(page_data + 4096 - offset - sizeof(slot.second), &slot.second, sizeof(slot.second));
            offset += sizeof(slot.second);
        }

        // Write the page_data buffer to the output stream
        out.write(page_data, sizeof(page_data));
    }

    // Function to read a page from a binary input stream
    bool read_from_data_file(istream &in) {
        char page_data[4096] = {0}; // Buffer to hold page data
        in.read(page_data, 4096); // Read data from input stream

        streamsize bytes_read = in.gcount();
        if (bytes_read == 4096) {
            //DONE
            //Process data to fill the records, slot_directory, and overflowPointerIndex
            //records
            int num, offset = 0;

            memcpy(&num, page_data + 4096 - sizeof(num) - offset, sizeof(num));
            offset += sizeof(num);

            memcpy(&overflowPointerIndex, page_data + 4096 - sizeof(overflowPointerIndex) - offset, sizeof(overflowPointerIndex));
            offset += sizeof(overflowPointerIndex);

            // slot_directory
            for (int i = 0; i < num; i++) {
                int first, second;

                memcpy(&first, page_data + 4096 - offset - sizeof(first), sizeof(first));
                offset += sizeof(first);

                memcpy(&second, page_data + 4096 - offset - sizeof(second), sizeof(second));
                offset += sizeof(second);
            
                cur_size =+ sizeof(int) * 2;
                slot_directory.emplace_back(first, second);
            }

            // processing records
            for (const auto& slot: slot_directory){
                offset = slot.first;
            
                int id, manager_id;
                memcpy(&id, page_data + offset, sizeof(int));
                offset += sizeof(int);
                memcpy(&manager_id, page_data + offset, sizeof(int));
                offset += sizeof(int);


                int name_len, bio_len;
                memcpy(&name_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                string name(page_data + offset, name_len);
                offset += name_len;
                memcpy(&bio_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                string bio(page_data + offset, bio_len);
                offset += bio_len;

                Record r(id, name, bio, manager_id);

                records.push_back(r);
                cur_size += r.get_size();
                records_size += r.get_size();
            }
            return true;
        }

        if (bytes_read > 0) {
            cerr << "Incomplete read: Expected 4096 bytes, but only read " << bytes_read << " bytes." << endl;
        }

        return false;
    }
};

class HashIndex {
private:
    const size_t maxCacheSize = 1; // Maximum number of pages in the buffer
    const int Page_SIZE = 4096; // Size of each page in bytes
    vector<int> PageDirectory{vector<int>(256,-1)}; // Map h(id) to a bucket location in EmployeeIndex(e.g., the jth bucket)
    // can scan to correct bucket using j*Page_SIZE as offset (using seek function)
    // can initialize to a size of 256 (assume that we will never have more than 256 regular (i.e., non-overflow) buckets)
    int nextFreePage; // Next place to write a bucket
    string fileName;

    // Function to compute hash value for a given ID
    int compute_hash_value(int id) {
        int hash_value;

        // DONE: Implement the hash function h = id mod 2^8
        hash_value = id % (1 << 8);

        return hash_value;
    }

    // Function to add a new record to an existing page in the index file
    void addRecordToIndex(int pageIndex,/* Page &page,*/ Record &record) {
        // Open index file in binary mode for updating
        fstream indexFile(fileName, ios::binary | ios::in | ios::out);

        if (!indexFile) {
            cerr << "Error: Unable to open index file for adding record." << endl;
            return;
        }

        Page page;
        // Seek to the appropriate position in the index file
        indexFile.seekp(pageIndex * Page_SIZE, ios::beg); 

        page.read_from_data_file(indexFile);

        // Check if the page has overflow
        while (page.overflowPointerIndex != -1) {
            pageIndex = page.overflowPointerIndex;
            indexFile.seekp(pageIndex * Page_SIZE, ios::beg);
            Page overflowPage;
            overflowPage.read_from_data_file(indexFile);
            page = std::move(overflowPage);
        }

        // DONE: Insert record to page and write data to file
        // Insert record to page
        if( page.insert_record_into_page(record)){
            indexFile.seekp(pageIndex * Page_SIZE, ios::beg);
            page.write_into_data_file(indexFile);
        } else {
            // If the record couldn't be inserted into the current page, insert it into the overflow page
            page.overflowPointerIndex = nextFreePage;
            nextFreePage++;
            Page overflowPage;
            if (!overflowPage.insert_record_into_page(record)) {
                cerr << "ERROR: couldn't insert into overflow page" << endl;
                return;
            }

            // Write page with its new overflow 
            indexFile.seekp(pageIndex * Page_SIZE, std::ios::beg);
            page.write_into_data_file(indexFile);

            // Write the overflow page
            indexFile.seekp(page.overflowPointerIndex * Page_SIZE, std::ios::beg);
            overflowPage.write_into_data_file(indexFile);
        }

        // Close the index file
        indexFile.close();
    }

    // Function to search for a record by ID in a given page of the index file
    void searchRecordByIdInPage(int pageIndex, int id) {
        // Open index file in binary mode for reading
        fstream indexFile(fileName, ios::binary | ios::out | ios::in);

        // Seek to the appropriate position in the index file
        indexFile.seekp(pageIndex * Page_SIZE, ios::beg);

        // Read the page from the index file
        Page page;
        page.read_from_data_file(indexFile);

        // DONE:
        //  - Search for the record by ID in the page
        for (auto& record: page.records) {
            if(record.id == id){
                record.print();
                return;
            }
        }
        //  - Check for overflow pages and report if record with given ID is not found
        if(page.overflowPointerIndex != -1){
            searchRecordByIdInPage(page.overflowPointerIndex, id);
            return;
        }
        cout << "Record with ID " << id << " not found." << endl;

        // Free the memory occupied by the page
        page.records.clear();
        page.slot_directory.clear();

        // Close the index file
        indexFile.close();

    }

public:
    HashIndex(string indexFileName) : nextFreePage(0), fileName(indexFileName) { }

    // Function to create hash index from Employee CSV file
    void createFromFile(string csvFileName) {
        // Read CSV file and add records to index
        // Open the CSV file for reading
        ifstream csvFile(csvFileName);

        // Start blank page     
        Page p;
        fstream indexFile;
        indexFile.open(fileName, ios::binary | ios::out | ios::in | ios::trunc);
        indexFile.seekp(0, ios::beg);
        p.write_into_data_file(indexFile);
        // Close the index file
        indexFile.close();

        string line;

        // Read each line from the CSV file
        while (getline(csvFile, line)) {
            // Parse the line and create a Record object
            stringstream ss(line);
            string item;
            vector<string> fields;
            while (getline(ss, item, ',')) {
                fields.push_back(item);
            }
            Record record(fields);

            // DONE:
            //   - Compute hash value for the record's ID using compute_hash_value() function.
            int hash_value = compute_hash_value(record.id);

            ///   - Get the page index from PageDirectory. If it's not in PageDirectory, define a new page using nextFreePage.
            int page_index = PageDirectory[hash_value];
            
            if(page_index == -1){
                //Create a new page
                Page p;
                page_index = nextFreePage;
                PageDirectory[hash_value] = page_index;
                nextFreePage++;
                p.insert_record_into_page(record);

                // write to spot in the file
                fstream indexFile(fileName, ios::binary | ios::out | ios::in);
                indexFile.seekp(page_index * Page_SIZE, ios::beg);
                p.write_into_data_file(indexFile);
                indexFile.close();
            } else {
                //   - Insert the record into the appropriate page in the index file using addRecordToIndex() function.
                addRecordToIndex(page_index, record);
            }
        }
        // Close the CSV file
        csvFile.close();
    }

    // Function to search for a record by ID in the hash index
    void findAndPrintEmployee(int id) {
        // Open index file in binary mode for reading
        ifstream indexFile(fileName, ios::binary | ios::in);

        //TEST
        if (!indexFile) {
            cerr << "Error: Unable to open index file '" << fileName << "' for searching." << endl;
            return;
        }

        // DONE:
        //  - Compute hash value for the given ID using compute_hash_value() function
        int hashValue = compute_hash_value(id);

        // Get the page index from PageDirectory
        int pageIndex;
        if (hashValue >= PageDirectory.size()) {
            // Create a new page
            PageDirectory.push_back(nextFreePage);
            pageIndex = nextFreePage;
            nextFreePage++;
            // cout << "Record with ID " << id << " not found." << endl;
            // return;
        } else {
            pageIndex = PageDirectory[hashValue];
        }

        //  - Search for the record in the page corresponding to the hash value using searchRecordByIdInPage() function
        searchRecordByIdInPage(pageIndex, id);

        // Close the index file
        indexFile.close();
    }
};


