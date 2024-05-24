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

    // // Constructor accepting an initializer_list<string>
    // Record(initializer_list<string> fields) {
    //     vector<string> field_vec(fields);
    //     id = stoi(field_vec[0]);
    //     name = field_vec[1];
    //     bio = field_vec[2];
    //     manager_id = stoi(field_vec[3]);
    // }

    // Constructor accepting a vector<string>
    Record(vector<string> &fields) {
        id = stoi(fields[0]);
        name = fields[1];
        bio = fields[2];
        manager_id = stoi(fields[3]);
    }

    // Function to get the size of the record
    int get_size() {
        return int(sizeof(int) * 2 + bio.length() + name.length());
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

    void print() const {                                //added const to call record.print()
        cout << "\tID: " << id << "\n";
        cout << "\tNAME: " << name << "\n";
        cout << "\tBIO: " << bio << "\n";
        cout << "\tMANAGER_ID: " << manager_id << "\n";
    }
};

class Page {
public:
    vector<Record> records; // Data_Area containing the records
    vector<pair<int, int>> slot_directory; // Slot directory containing offset and size of each record
    int cur_size = 0; // Current size of the page
    int overflowPointerIndex; // Offset of overflow page, set to -1 by default

    // Constructor
    Page() : overflowPointerIndex(-1) {}

    // Function to insert a record into the page
    bool insert_record_into_page(Record r) {
        if (cur_size + r.get_size() >= 4096) { // Check if page size limit exceeded
            return false; // Cannot insert the record into this page
        } else {
            records.push_back(r);
            cur_size += r.get_size();
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

        // DONE:
        //  - Write slot_directory in reverse order into page_data buffer.
        for (int i = slot_directory.size() - 1; i >= 0; i--) {
            const auto &entry = slot_directory[i];
            memcpy(page_data + offset, &entry.first, sizeof(int));
            offset += sizeof(int);
            memcpy(page_data + offset, &entry.second, sizeof(int));
            offset += sizeof(int);
        }
        //  - Write overflowPointerIndex into page_data buffer.
        memcpy(page_data + offset, &overflowPointerIndex, sizeof(int));
        offset += sizeof(int);
        //  You should write the first entry of the slot_directory, which have the info about the first record at the bottom of the page, before overflowPointerIndex.

        // Write the page_data buffer to the output stream
        out.write(page_data, sizeof(page_data));
    }

    // Function to read a page from a binary input stream
    bool read_from_data_file(istream &in) {
        char page_data[4096] = {0}; // Buffer to hold page data
        in.read(page_data, 4096); // Read data from input stream

        streamsize bytes_read = in.gcount();
        if (bytes_read == 4096) {
            // DONE: Process data to fill the records, slot_directory, and overflowPointerIndex
            int offset = 0;
            int slot_offset, slot_size;

           // Read records from page_data buffer
            while (offset < 4096 - sizeof(int)) {
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

                // Create a vector<string> and a Record object using the vector constructor
                vector<string> fields;
                fields.assign({to_string(id), name, bio, to_string(manager_id)});
                records.emplace_back(Record(fields));
            }
            // Read slot_directory from page_data buffer
            while (offset < 4096 - sizeof(int)) {
                slot_offset = *reinterpret_cast<int*>(page_data + offset);
                offset += sizeof(int);
                slot_size = *reinterpret_cast<int*>(page_data + offset);
                offset += sizeof(int);
                slot_directory.emplace_back(slot_offset, slot_size);
            }

            // Read overflowPointerIndex from page_data buffer
            overflowPointerIndex = *reinterpret_cast<int*>(page_data + offset);

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
    vector<int> PageDirectory; // Map h(id) to a bucket location in EmployeeIndex(e.g., the jth bucket)
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
    void addRecordToIndex(int pageIndex, Page &page, Record &record) {
        // Open index file in binary mode for updating
        fstream indexFile(fileName, ios::binary | ios::in | ios::out);

        if (!indexFile) {
            cerr << "Error: Unable to open index file for adding record." << endl;
            return;
        }

        //TEST
        //cout << "Adding record with ID: " << record.id << endl;

        // Check if the page has overflow
        if (page.overflowPointerIndex == -1) {
            // DONE: Create overflow page using nextFreePage. update nextFreePage index and pageIndex
            //TEST
            //cout << "Creating overflow page at index: " << nextFreePage << endl;
            page.overflowPointerIndex = nextFreePage;
            PageDirectory.push_back(nextFreePage);
            nextFreePage++;
        }

        // Seek to the appropriate position in the index file
        indexFile.seekp(pageIndex * Page_SIZE, ios::beg);

        // DONE: Insert record to page and write data to file
        // Insert record to page
        if (!page.insert_record_into_page(record)) {
            //TEST 
            //cout << "Record could not be inserted into main page, inserting into overflow page." << endl;
            // If the record couldn't be inserted into the current page, insert it into the overflow page
            int overflowPageIndex = page.overflowPointerIndex;
            Page overflowPage;
            indexFile.seekg(overflowPageIndex * Page_SIZE, ios::beg);
            overflowPage.read_from_data_file(indexFile);
            overflowPage.insert_record_into_page(record);

            // Write the overflow page to the file
            indexFile.seekp(overflowPageIndex * Page_SIZE, ios::beg);
            overflowPage.write_into_data_file(indexFile);

            // Free the memory occupied by the overflow page
            overflowPage.records.clear();
            overflowPage.slot_directory.clear();
        }
        //else {
            //TEST
            //cout << "Record inserted into main page." << endl;
        //}
        // Write the page to the file
        indexFile.seekp(pageIndex * Page_SIZE, ios::beg);
        page.write_into_data_file(indexFile);

        // Free the memory occupied by the page
        page.records.clear();
        page.slot_directory.clear();

        // Close the index file
        indexFile.close();
    }

    // Function to search for a record by ID in a given page of the index file
    void searchRecordByIdInPage(int pageIndex, int id) {
        // Open index file in binary mode for reading
        ifstream indexFile(fileName, ios::binary | ios::in);

        //TEST
        if (!indexFile) {
            cerr << "Error: Unable to open index file '" << fileName << "' for searching." << endl;
            return;
        }

        // Seek to the appropriate position in the index file
        indexFile.seekg(pageIndex * Page_SIZE, ios::beg);

        // Read the page from the index file
        Page page;
        page.read_from_data_file(indexFile);

        // DONE:
        //  - Search for the record by ID in the page
        bool found = false;
        for (const auto& record : page.records) {
            if (record.id == id) {
                record.print();
                found = true;
                break;
            }
        }
        //  - Check for overflow pages and report if record with given ID is not found
        if (!found && page.overflowPointerIndex != -1) {
            int overflowPageIndex = page.overflowPointerIndex;
            searchRecordByIdInPage(overflowPageIndex, id);
        } else if (!found) {
            cout << "Record with ID " << id << " not found." << endl;
        }

        // Free the memory occupied by the page
        page.records.clear();
        page.slot_directory.clear();

        // Close the index file
        indexFile.close();
    }

public:
    HashIndex(string indexFileName) : nextFreePage(0), fileName(indexFileName) {
    }

    // Function to create hash index from Employee CSV file
    void createFromFile(string csvFileName) {
        // Read CSV file and add records to index
        // Open the CSV file for reading
        ifstream csvFile(csvFileName);

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

            // TEST: Print the details of the Record object
            // cout << "Record details:" << endl;
            // record.print();
            // cout << endl;

            // DONE:
            //   - Compute hash value for the record's ID using compute_hash_value() function.
            int hashValue = compute_hash_value(record.id);

            //   - Get the page index from PageDirectory. If it's not in PageDirectory, define a new page using nextFreePage.
            int pageIndex;
            if (hashValue >= PageDirectory.size()) {
                // Create a new page
                PageDirectory.push_back(nextFreePage);
                pageIndex = nextFreePage;
                nextFreePage++;
            } else {
                pageIndex = PageDirectory[hashValue];
            }

            // Check if the index file exists
            bool fileExists = std::ifstream(fileName).good();

            // Open the index file in binary mode for reading and writing
            fstream indexFile;
            if (fileExists){
                indexFile.open(fileName, ios::binary | ios::in | ios::out);
            } else {
                indexFile.open(fileName, ios::binary | ios::out); // Create a new file
            }

            // Create a Page object
            Page page;

            // Read the page from the index file (if it exists)
            if (fileExists) {
                indexFile.seekg(pageIndex * Page_SIZE, ios::beg);
                page.read_from_data_file(indexFile);
            }

            //   - Insert the record into the appropriate page in the index file using addRecordToIndex() function.
            addRecordToIndex(pageIndex, page, record);

            // Close the index file
            indexFile.close();

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
