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

    Record(vector<string> &fields) {
        id = stoi(fields[0]);
        name = fields[1];
        bio = fields[2];
        manager_id = stoi(fields[3]);
    }

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

class page{ // Take a look at Figure 9.7 and read Section 9.6.2 [Page Organization for Variable Length Records] 
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
            records.push_back(r);
            slot_directory.push_back(make_pair(cur_size, record_size)); // added
            cur_size += r.get_size(); 
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
        if (bytes_read != 4096) {
            if (bytes_read > 0) {
                cerr << "Incomplete read: Expected " << 4096 << " bytes, but only read " << bytes_read << " bytes." << endl;
            }
            return false;
        }
        else {
            int offset = 0;
            // Process data to fill the slot directory and the records to handle it according to the structure
            // Assuming slot directory is processed here or elsewhere depending on your serialization method
            while (offset < 4096) {
                int id, manager_id;
                string name, bio;

                // Read the fixed-size fields.
                memcpy(&id, page_data + offset, sizeof(int));
                offset += sizeof(int);
                memcpy(&manager_id, page_data + offset, sizeof(int));
                offset += sizeof(int);

                // Read the variable-length fields.
                int name_len, bio_len;
                memcpy(&name_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                name.assign(page_data + offset, name_len);
                offset += name_len;
                memcpy(&bio_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                bio.assign(page_data + offset, bio_len);
                offset += bio_len;

                // Create a new record and add it to the page.
                vector<string> fields;
                fields.push_back(to_string(id));
                fields.push_back(name);
                fields.push_back(bio);
                fields.push_back(to_string(manager_id));
                records.push_back(Record(fields));
                slot_directory.push_back(make_pair(offset - name_len - bio_len - 2 * sizeof(int), name_len + bio_len + 2 * sizeof(int)));
            }
                        
            return true;
        }
    }
};

class StorageManager {

public:
    string filename;  // Name of the file (EmployeeRelation.dat) where we will store the Pages 
    fstream data_file; // fstream to handle both input and output binary file operations
    vector <page> buffer; // You can have maximum of 3 Pages.
    
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

        // addeding check 
        if (!csvFile.is_open()) {
            cerr << "Failed to open CSV file: " << csvFilename << endl;
            return;  // Return if the CSV file cannot be opened
        }
        
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
            if (page_number < buffer.size() && !buffer[page_number].insert_record_into_page(r)) {

                // If the current page is full, move to the next page
                page_number++;
 
                // If you have used all your 3 pages and they are filled up with records you may use write_into_data_file() to write the pages into the data file
                // It is like dumping all the information in you pages (Which are in your main memory) into the .dat file (which is in External Storage)
                if (page_number >= buffer.size()) {
                    for (page& p : buffer) {
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
        csvFile.close();  // Close the CSV file
    }

    // Searches for an Employee by ID in the binary data_file using the page and prints if found
    void findAndPrintEmployee(int searchId) {
        //adding another check 
        if (!data_file.is_open()) {
            cerr << "Data file is not open." << endl;
        return;
    }
        data_file.seekg(0, ios::beg);  // Rewind the data_file to the beginning for reading

        /*** TO_DO ***/
        // Use [ read_from_data_file(data_file)] to read a page from your datafile to the page buffer[3];
        // You can read into these 3 pages at once or use the following loop to read pages one by one. 
        buffer.resize(3); // You can have a maximum of 3 Pages
        int page_number = 0;
        bool found = false;
        while(buffer[page_number].read_from_data_file(data_file)){

            // Now that you have a page loaded with the data you stored previously, use the slot directory to find the desired id.
            // If you have not found the id in these 1/3 pages, repeat the process: reset the pages[3] and load next 3 pages from the EmployeeRelation.dat file   
            // Search for the record with the given ID in the current page
            for (const Record& r : buffer[page_number].records) {
                if (r.id == searchId) {
                    r.print();
                    found = true;
                    break;
                }
            }

            // Move to the next page
            page_number++;
            if (page_number >= buffer.size()) {
                page_number = 0;
                buffer.clear();
                buffer.resize(3);
            }
        }
      
        // Print not found message if no match from any of the records
        if (!found) {
            cout << "Employee with ID " << searchId << " not found." << endl;
        }
    }
};
