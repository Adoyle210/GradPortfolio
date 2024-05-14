/*** This is just a Skeleton/Starter Code for the External Storage Assignment. This is by no means absolute, in terms of assignment approach/ used functions, etc. ***/
/*** You may modify any part of the code, as long as you stick to the assignments requirements we do not have any issue ***/

// Include necessary standard library headers
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <bitset>
#include <fstream> // Include this header
using namespace std; // Include the standard namespace

class Record {
public:
    int id, manager_id; // Employee ID and their manager's ID
    std::string bio, name; // Fixed length string to store employee name and biography

    Record(vector<std::string> &fields) {
        id = std::stoi(fields[0]);
        name = fields[1];
        bio = fields[2];
        manager_id = std::stoi(fields[3]);
    }

    //added 
    Record(int id, std::string name, std::string bio, int manager_id)
        : id(id), name(std::move(name)), bio(std::move(bio)), manager_id(manager_id) {}


    //You may use this for debugging / showing the record to standard output. 
    void print() const {
        cout << "\tID: " << id << "\n";
        cout << "\tNAME: " << name << "\n";
        cout << "\tBIO: " << bio << "\n";
        cout << "\tMANAGER_ID: " << manager_id << "\n";
    }

    int get_size(){
        return int(sizeof(int) * 2 + bio.length() + name.length());
    }
     // Serialize the record for writing to file
    string serialize() const {
        ostringstream oss;
        oss.write(reinterpret_cast<const char*>(&id), sizeof(id));
        oss.write(reinterpret_cast<const char*>(&manager_id), sizeof(manager_id));
        int name_len = name.size();
        int bio_len = bio.size();
        oss.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
        oss.write(name.c_str(), name.size());
        oss.write(reinterpret_cast<const char*>(&bio_len), sizeof(bio_len));
        oss.write(bio.c_str(), bio.size());
        return oss.str();
    }
};

class Page{ // Take a look at Figure 9.7 and read Section 9.6.2 [Page Organization for Variable Length Records] 
public:
    vector <Record> records; // This is the Data_Area, which contains the records. 
    vector <pair <int, int> > slot_directory; // This slot directory contains the starting position (offset), and size of the record. 
                                             // The starting position i refers to the location of record i in the records and size of that particular record i. 
                                             // This information is crucial, you can load record using this if you know if your record is starting from, e.g., 
                                            // byte 128 and has size 70, you read all characters from [128 to 128 + 70]
    
    int cur_size = 0; // holds the current size of the 

    // Write a Record into your page
    bool insert_record_into_page(Record r){
        int record_size = r.get_size();
        // Take a look at Figure 9.9 and read the Section 9.7.2 [Record Organization for Variable Length Records]
        // You may adopt any of the approaches mentioned their. E.g., id $ name $ bio $ manager_id $ separating records with a delimiter / the alternative approaches
        if(cur_size + record_size >= 4096){     // Check the current size of your page. Your page has 4KB memory for storing the records and the slot directory information.   

            //You cannot insert the current record into this page
            return false;
        }
        else{
            int offset = cur_size;
            records.push_back(r);
            slot_directory.push_back(make_pair(offset, record_size));
            cur_size += record_size;
            // update slot directory information
            // slot_directory.push_back(make_pair(offset, serialized.size()));

            return true;
        }
        

    }

    // Function to write the page to a binary output stream, i.e., EmployeeRelation.dat file
    void write_into_data_file(ostream& out) const {
        //Write the records and slot directory information into your data file. You are basically writing 4KB into the datafile. 
        //You must maintain a fixed size of 4KB so there may be some unused empty spaces. 
        
        char page_data[4096] = {0}; // Let's write all the information of the page into this char array. So that we can write the page into the data file in one go.
        int offset = 0;

        // If you look at figure 9.7, you'll find that there are spaces allocated for the slot-directory. 
        // You can structure your page in your own way, such as allocate first x bytes of memory to store the slot-directory information
        //  sizeof(int) bytes to parse these numbers. 
        // After those x bytes, you start storing your records. 
        // You can definitely use $ (delimiter) while storing the slot directory informations /(or,) as you know that these are integers(sizeof(int)) you can read 
        

        for (const auto& record : records) {
            string serialized = record.serialize();

            memcpy(page_data + offset, serialized.c_str(), serialized.size());

            offset += serialized.size();
        }
        //the above loop just read the id, name, bio etc. You'll also need to store the slot-directory information. So that you can use the slot-directory
        // to retrieve a record.
        for (const auto& slots : slot_directory) {
            // insert the slot directory information into the page_data
            memcpy(page_data + offset, &slots.first, sizeof(int)); // Write the offset.
            offset += sizeof(int);
            memcpy(page_data + offset, &slots.second, sizeof(int)); // Write the size.
            offset += sizeof(int);
        }

        // Write the page_data to the EmployeeRelation.dat file 
        out.write(page_data, sizeof(page_data)); // Always write exactly 4KB
    }

    // Read a page from a binary input stream, i.e., EmployeeRelation.dat file to populate a page object
    bool read_from_data_file(istream& in) {
        // Read all the records and slot_directory information from your .dat file
        // Remember you are reading a chunk of 4098 byte / 4 KB from the data file to your main memory. 
        char page_data[4096] = {0};
        in.read(page_data, 4096);

        streamsize bytes_read = in.gcount();
        // You may populate the records and slot_directory from the 4 KB data you just read.
        //addomg for testing
       if (bytes_read == 4096) {
            int offset = 0;

            // Process data to fill the slot directory and the records to handle it according to the structure
            // Assuming slot directory is processed here or elsewhere depending on your serialization method
            int slotDirectorySize;
            memcpy(&slotDirectorySize, page_data, sizeof(int));
            offset += sizeof(int);
            slot_directory.resize(slotDirectorySize / sizeof(pair<int, int>));
            memcpy(slot_directory.data(), page_data + offset, slotDirectorySize);

            offset += slotDirectorySize;
            // Read records
            for (const auto& [recordOffset, recordSize] : slot_directory) {
                string serialized(page_data + offset, recordSize);
                int id, managerId, nameLen, bioLen;
                memcpy(&id, serialized.c_str(), sizeof(int));
                memcpy(&managerId, serialized.c_str() + sizeof(int), sizeof(int));
                memcpy(&nameLen, serialized.c_str() + sizeof(int) * 2, sizeof(int));
                std::string name(serialized.c_str() + sizeof(int) * 3, nameLen);
                memcpy(&bioLen, serialized.c_str() + sizeof(int) * 3 + nameLen, sizeof(int));
                std::string bio(serialized.c_str() + sizeof(int) * 4 + nameLen, bioLen);
                records.push_back(Record(id, std::move(name), std::move(bio), managerId));
                offset += recordSize;
            }

            cur_size = offset;            
            return true;
        }
        if (bytes_read > 0) {
            cerr << "Incomplete read: Expected " << sizeof(page_data) << " bytes, but only read " << bytes_read << " bytes." << endl;
        }

        return false;
    }
};

class StorageManager {
public:
    string filename;  // Name of the file (EmployeeRelation.dat) where we will store the Pages 
    fstream data_file; // fstream to handle both input and output binary file operations
    vector <Page> buffer; // You can have maximum of 3 Pages.
    
    // Constructor that opens a data file for binary input/output; truncates any existing data file
    StorageManager(const string& filename) : filename(filename) {
        data_file.open(filename, ios::binary | ios::out | ios::in | ios::trunc);
        if (!data_file.is_open()) {  // Check if the data_file was successfully opened
            cerr << "Failed to open data_file: " << filename << endl;
            exit(EXIT_FAILURE);  // Exit if the data_file cannot be opened
        }
    }

    // Destructor closes the data file if it is still open
    ~StorageManager() {
        if (data_file.is_open()) {
            data_file.close();
        }
    }

    // Reads data from a CSV file and writes it to a binary data file as Employee objects
    void createFromFile(const string& csvFilename) {
        buffer.resize(3); // You can have maximum of 3 Pages.

        ifstream csvFile(csvFilename);  // Open the Employee.csv file for reading
        
        string line, name, bio;
        int id, manager_id;
        int page_number = 0; // Current page we are working on [at most 3 pages]

        // Read each line from the CSV file, parse it, and create Employee objects
        while (getline(csvFile, line)) {
            stringstream ss(line);
            string item;
            vector<string> fields;
            // Split line by commas
            while (getline(ss, item, ',')) {
                fields.push_back(item);
            }
            //create A record object            
            Record r = Record(fields);
            // r.print();
            /***TO_DO***/ 

            // Now we will insert that record information to page, i.e., buffer[page_number].
            // If the page has < 4KB data
                // We Insert the record into current page in the buffer
            // else:
                // go to the next page [page_number++]
                // Remember You only have 3 pages.
            //edit check to now have: page_number < buffer.size() 
            if (page_number < buffer.size()) {
                if (!buffer[page_number].insert_record_into_page(r)) {

                    // If the current page is full, move to the next page
                    page_number++;
    
                    // If you have used all your 3 pages and they are filled up with records you may use write_into_data_file() to write the pages into the data file
                    // It is like dumping all the information in you pages (Which are in your main memory) into the .dat file (which is in External Storage)
                    if (page_number >= buffer.size()) {
                        for (Page& p : buffer) {
                            p.write_into_data_file(data_file);
                        }
                        // Reset/Free the pages and start filling them up with the future records
                        page_number = 0;  
                        // page.clear();
                        buffer.clear();
                        buffer.resize(3);
                    }
                    // Insert the record into the new current page
                    buffer[page_number].insert_record_into_page(r);
                }
            }
            
        }
        csvFile.close();  // Close the CSV file
    }

    // Searches for an Employee by ID in the binary data_file using the page and prints if found
    void findAndPrintEmployee(int searchId) {
        data_file.seekg(0, ios::beg);  // Rewind the data_file to the beginning for reading
        buffer.resize(3);
        /*** TO_DO ***/
        // Use [ read_from_data_file(data_file)] to read a page from your datafile to the page buffer[3];
        // You can read into these 3 pages at once or use the following loop to read pages one by one. 
        int page_number = 0;

        while (true) {
            if (page_number < buffer.size()) {
                if (!buffer[page_number].read_from_data_file(data_file)) {
                    break;
                }
            } else {
                buffer[0].read_from_data_file(data_file);
                page_number = 1;
            }

            for (const Page& page : buffer) {
                for (const auto& [offset, size] : page.slot_directory) {
                    const std::string& serialized = page.records[offset].serialize();
                    int recordId;
                    memcpy(&recordId, serialized.c_str(), sizeof(int));
                    if (recordId == searchId) {
                        page.records[offset].print();
                        return;
                    }
                }
            }

            page_number++;
        }
      
        // Print not found message if no match from any of the records
        cout << "Employee with ID " << searchId << " not found." << endl;
    }
};