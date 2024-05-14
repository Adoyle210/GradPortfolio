/*** This is just a Skeleton/Starter Code for the External Storage Assignment. This is by no means absolute, in terms of assignment approach/ used functions, etc. ***/
/*** You may modify any part of the code, as long as you stick to the assignments requirements we do not have any issue ***/

// Include necessary standard library headers
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

class Record {
public:
    int id, manager_id;// Employee ID and their manager's ID
    string bio, name;// Fixed length string to store employee name and biography

    Record(vector<string>& fields) {
        try {
            id = stoi(fields[0]);
            name = fields[1];
            bio = fields[2];
            manager_id = stoi(fields[3]);
        } catch (const invalid_argument& e) {
            cerr << "Error parsing record: " << e.what() << endl;
        }
    }

    Record(int id, string name, string bio, int manager_id)
        : id(id), name(move(name)), bio(move(bio)), manager_id(manager_id) {}

    //You may use this for debugging / showing the record to standard output. 
    void print() const {
        cout << "\tID: " << id << "\n";
        cout << "\tNAME: " << name << "\n";
        cout << "\tBIO: " << bio << "\n";
        cout << "\tMANAGER_ID: " << manager_id << "\n";
    }

    int get_size() const {
        return sizeof(int) * 2 + bio.length() + name.length();
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

class Page {// Take a look at Figure 9.7 and read Section 9.6.2 [Page Organization for Variable Length Records] 
public:
    vector<Record> records; // This is the Data_Area, which contains the records. 
    vector<pair<int, int> > slot_directory;// This slot directory contains the starting position (offset), and size of the record. 
                                             // The starting position i refers to the location of record i in the records and size of that particular record i. 
                                             // This information is crucial, you can load record using this if you know if your record is starting from, e.g., 
                                            // byte 128 and has size 70, you read all characters from [128 to 128 + 70]
    
    int curSize = 0;// holds the current size of the 

    bool insert_record_into_page(const Record& r) {
        int recordSize = r.get_size();
        // Take a look at Figure 9.9 and read the Section 9.7.2 [Record Organization for Variable Length Records]
        // You may adopt any of the approaches mentioned their. E.g., id $ name $ bio $ manager_id $ separating records with a delimiter / the alternative approaches
        if (curSize + recordSize >= 4096) { // Check the current size of your page. Your page has 4KB memory for storing the records and the slot directory information.   
            //You cannot insert the current record into this page
            return false;
        } else {
            int offset = curSize;
            records.push_back(r);
            slot_directory.push_back(make_pair(offset, recordSize));
            curSize += recordSize;
            // update slot directory information
            // slot_directory.push_back(make_pair(offset, serialized.size()));
            return true;
        }
    }

    // Function to write the page to a binary output stream, i.e., EmployeeRelation.dat file
    void write_into_data_file(ostream& out) const {
        //Write the records and slot directory information into your data file. You are basically writing 4KB into the datafile. 
        //You must maintain a fixed size of 4KB so there may be some unused empty spaces. 
        
        char page_data[4096] = {0};// Let's write all the information of the page into this char array. So that we can write the page into the data file in one go.
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

        int slot_directory_size = slot_directory.size() * sizeof(pair<int, int>);
        memcpy(page_data, &slot_directory_size, sizeof(int));
        offset = sizeof(int);
        memcpy(page_data + offset, slot_directory.data(), slot_directory_size);
        offset += slot_directory_size;

        out.write(page_data, sizeof(page_data));
    }

    // Read a page from a binary input stream, i.e., EmployeeRelation.dat file to populate a page object
    bool read_from_data_file(istream& in) {
        // Read all the records and slot_directory information from your .dat file
        // Remember you are reading a chunk of 4098 byte / 4 KB from the data file to your main memory. 
        char page_data[4096] = {0};
        in.read(page_data, sizeof(page_data));
        streamsize bytesRead = in.gcount();
        // You may populate the records and slot_directory from the 4 KB data you just read.
        //addomg for testing
        if (bytesRead == sizeof(page_data)) {
            int offset = 0;
            int slot_directory_size;
            memcpy(&slot_directory_size, page_data, sizeof(int));
            offset += sizeof(int);

            int numSlots = slot_directory_size / sizeof(pair<int, int>);
            slot_directory.resize(numSlots);
            memcpy(slot_directory.data(), page_data + offset, slot_directory_size);
            offset += slot_directory_size;

            for (const auto& [recordOffset, recordSize] : slot_directory) {
                const string serialized(page_data + offset, recordSize);
                int id, managerId, nameLen, bioLen;
                memcpy(&id, serialized.c_str(), sizeof(int));
                memcpy(&managerId, serialized.c_str() + sizeof(int), sizeof(int));
                memcpy(&nameLen, serialized.c_str() + sizeof(int) * 2, sizeof(int));
                string name(serialized.c_str() + sizeof(int) * 3, nameLen);
                memcpy(&bioLen, serialized.c_str() + sizeof(int) * 3 + nameLen, sizeof(int));
                string bio(serialized.c_str() + sizeof(int) * 4 + nameLen, bioLen);
                records.push_back(Record(id, move(name), move(bio), managerId));
                offset += recordSize;
            }

            curSize = offset;
            return true;
        }

        if (bytesRead > 0) {
            cerr << "Incomplete read: Expected " << sizeof(page_data) << " bytes, but only read " << bytesRead << " bytes." << endl;
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
    StorageManager(const string& filename)
        : filename(filename) {
        data_file.open(filename, ios::binary | ios::out | ios::in | ios::trunc);
        if (!data_file.is_open()) {
            cerr << "Failed to open data file: " << filename << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Destructor closes the data file if it is still open
    ~StorageManager() {
        if (data_file.is_open()) {
            data_file.close();
        }
    }

    // Reads data from a CSV file and writes it to a binary data file as Employee objects
    void create_from_file(const string& csvFilename) {
        buffer.resize(3); // You can have maximum of 3 Pages.
        ifstream csvFile(csvFilename); // Open the Employee.csv file for reading


        string line;
        int pageNumber = 0;

        // Read each line from the CSV file, parse it, and create Employee objects
        while (getline(csvFile, line)) {
            stringstream ss(line);
            string item;
            vector<string> fields;
            while (getline(ss, item, ',')) {
                fields.push_back(item);
            }
            Record r = Record(fields);

            if (pageNumber < buffer.size()) {
                if (!buffer[pageNumber].insert_record_into_page(r)) {
                    pageNumber++;

                    if (pageNumber >= buffer.size()) {
                        for (Page& p : buffer) {
                            p.write_into_data_file(data_file);
                        }
                        pageNumber = 0;
                        buffer.clear();
                        buffer.resize(3);
                    }
                    buffer[pageNumber].insert_record_into_page(r);
                }
            }
        }
        csvFile.close();

        for (Page& p : buffer) {
            p.write_into_data_file(data_file);
        }
    }
    // Searches for an Employee by ID in the binary data_file using the page and prints if found
    void findAndPrintEmployee(int searchId) {
        data_file.seekg(0, ios::beg);
        buffer.resize(3);
        int pageNumber = 0;

        while (true) {
            if (pageNumber < buffer.size()) {
                if (!buffer[pageNumber].read_from_data_file(data_file)) {
                    break;
                }
            } else {
                buffer[0].read_from_data_file(data_file);
                pageNumber = 1;
            }

            for (const Page& page : buffer) {
                for (const auto& [offset, size] : page.slot_directory) {
                    const string& serialized = page.records[offset].serialize();
                    if (!serialized.empty()) {
                        int recordId;
                        memcpy(&recordId, serialized.c_str(), sizeof(int));
                        if (recordId == searchId) {
                            page.records[offset].print();
                            return;
                        }
                    }
                }
            }

            pageNumber++;
        }
        
        // Print not found message if no match from any of the records
        cout << "Employee with ID " << searchId << " not found." << endl;
    }
};
