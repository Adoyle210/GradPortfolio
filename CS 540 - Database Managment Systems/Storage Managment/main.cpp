/*
Skeleton code for storage management
*/

#include <string>
#include <ios>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include "classes.h"
using namespace std;

int main(int argc, char* argv[]) {
    // Initialize the Storage Manager Class with the Binary .dat file name we want to create
    StorageManager manager("EmployeeRelation.dat");

    // Assuming the Employee.CSV file is in the same directory, 
    // we want to read from the Employee.csv and write into the new data_file
    manager.createFromFile("Employee.csv");

    // Searching for Employee IDs Using [manager.findAndPrintEmployee(id)]
    /***TO_DO***/ 

    //using argc
    // if (argc < 2) {
    //     cerr << "Error Usage: " << argv[0] << " <employee_id> " << endl;
    //     return 1;
    // }

    // int searchID = stoi(argv[1]);
    // manager.findAndPrintEmployee(searchID);

    //allowing multiple searches with a loop
    int searchID;
    bool keepSearching = true;

    while (keepSearching) {
        cout << "Enter an employee ID to search for (or 'q' to quit): ";
        string input;
        getline(cin, input);

        if (input == "q") {
            keepSearching = false;
            continue;
        }

        try {
            searchID = stoi(input);
            manager.findAndPrintEmployee(searchID);
        } catch (const invalid_argument& e) {
            cerr << "Invalid input. Please enter a valid employee ID or 'q' to quit." << endl;
        }catch (const out_of_range& e) {
            cout << "Employee with ID " << searchID << " not found." << endl;
        }
    }

    return 0;
}