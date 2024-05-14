#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class Record {
public:
    int id, manager_id;
    string bio, name;

    Record(const vector<string>& fields) {
        id = stoi(fields[0]);
        name = fields[1];
        bio = fields[2];
        manager_id = stoi(fields[3]);
    }

    void print() const {
        cout << "\tID: " << id << "\n";
        cout << "\tNAME: " << name << "\n";
        cout << "\tBIO: " << bio << "\n";
        cout << "\tMANAGER_ID: " << manager_id << "\n";
    }

    int get_size() const {
        return sizeof(int) * 2 + bio.length() + name.length();
    }

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

class Page {
public:
    vector<Record> records;
    vector<pair<int, int> > slot_directory;
    int cur_size = 0;

    bool insert_record_into_page(const Record& r) {
        int record_size = r.get_size();
        if (cur_size + record_size >= 4096) {
            return false;
        } else {
            records.push_back(r);
            slot_directory.push_back(make_pair(cur_size, record_size));
            cur_size += record_size;
            return true;
        }
    }

    void write_into_data_file(ostream& out) const {
        char page_data[4096] = {0};
        int offset = 0;

        for (const auto& record : records) {
            string serialized = record.serialize();
            memcpy(page_data + offset, serialized.c_str(), serialized.size());
            offset += serialized.size();
        }

        for (const auto& slots : slot_directory) {
            memcpy(page_data + offset, &slots.first, sizeof(int));
            offset += sizeof(int);
            memcpy(page_data + offset, &slots.second, sizeof(int));
            offset += sizeof(int);
        }

        out.write(page_data, sizeof(page_data));
    }

    bool read_from_data_file(istream& in) {
        char page_data[4096] = {0};
        in.read(page_data, 4096);

        streamsize bytes_read = in.gcount();
        if (bytes_read != 4096) {
            cerr << "Incomplete read: Expected 4096 bytes, but only read " << bytes_read << " bytes." << endl;
            return false;
        } else {
            int offset = 0;
            while (offset < 4096) {
                int id, manager_id;
                string name, bio;

                memcpy(&id, page_data + offset, sizeof(int));
                offset += sizeof(int);
                memcpy(&manager_id, page_data + offset, sizeof(int));
                offset += sizeof(int);

                int name_len, bio_len;
                memcpy(&name_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                name.assign(page_data + offset, name_len);
                offset += name_len;
                memcpy(&bio_len, page_data + offset, sizeof(int));
                offset += sizeof(int);
                bio.assign(page_data + offset, bio_len);
                offset += bio_len;

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
    string filename;
    fstream data_file;
    vector<Page> buffer;

    StorageManager(const string& filename) : filename(filename) {
        data_file.open(filename, ios::binary | ios::out | ios::in | ios::trunc);
        if (!data_file.is_open()) {
            cerr << "Failed to open data_file: " << filename << endl;
            exit(EXIT_FAILURE);
        }
    }

    ~StorageManager() {
        if (data_file.is_open()) {
            data_file.close();
        }
    }

    void createFromFile(const string& csvFilename) {
        buffer.resize(3);

        ifstream csvFile(csvFilename);
        if (!csvFile.is_open()) {
            cerr << "Failed to open CSV file: " << csvFilename << endl;
            return;
        }
        
        string line;
        while (getline(csvFile, line)) {
            stringstream ss(line);
            string item;
            vector<string> fields;
            while (getline(ss, item, ',')) {
                fields.push_back(item);
            }
            Record r(fields);

            int page_number = 0;
            while (page_number < buffer.size() && !buffer[page_number].insert_record_into_page(r)) {
                page_number++;
            }

            if (page_number >= buffer.size()) {
                for (Page& p : buffer) {
                    p.write_into_data_file(data_file);
                }
                page_number = 0;  
                buffer.clear();
                buffer.resize(3);
            }
            buffer[page_number].insert_record_into_page(r);
        }
        csvFile.close();
        // After finishing writing, close the data file
        data_file.close();
    }

    void findAndPrintEmployee(int searchId) {
        // Reopen the data file before reading from it
        data_file.open(filename, ios::binary | ios::in);
        if (!data_file.is_open()) {
            cerr << "Failed to open data_file: " << filename << endl;
            exit(EXIT_FAILURE);
        }

        buffer.resize(3);
        int page_number = 0;
        bool found = false;
        while(buffer[page_number].read_from_data_file(data_file)) {
            for (const Record& r : buffer[page_number].records) {
                if (r.id == searchId) {
                    r.print();
                    found = true;
                    break;
                }
            }
            page_number++;
            if (page_number >= buffer.size()) {
                page_number = 0;
                buffer.clear();
                buffer.resize(3);
            }
        }
      
        if (!found) {
            cout << "Employee with ID " << searchId << " not found." << endl;
        }
        // After finishing reading, close the data file
        data_file.close();
    }
};

int main() {
    StorageManager manager("EmployeeRelation.dat");
    manager.createFromFile("Employee.csv");
    manager.findAndPrintEmployee(11432121);
    return 0;
}
