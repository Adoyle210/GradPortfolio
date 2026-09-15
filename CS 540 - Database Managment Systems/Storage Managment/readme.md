# Storage Management

Assignment for storing and searching an Employee relation in an external data file.

## Summary
This project reads `Employee.csv`, stores the records in a binary relation file, and supports searching by employee id.

## Tech Stack
- C++
- Binary file I/O
- CSV parsing

## Assignment
Assume a relation `Employee(id, name, bio, manager-id)`. The values of `id` and `manager-id` are integers with fixed sizes of 8 bytes. The values of `name` and `bio` are character strings with maximum lengths of 200 and 500 characters, respectively.

### Input File
The program reads the input Employee relation from `Employee.csv` in the current working directory.

### Data File Creation
The program stores the records in a new binary data file named `EmployeeRelation.dat` in the current working directory.

### Searching the Data File
After file creation, the program accepts an Employee id on the command line and searches `EmployeeRelation.dat` for matching records.

## Build and Run
Compile with:
```sh
g++ -std=c++11 main.cpp -o main.out
```

Run with:
```sh
./main.out
```

## Notes
Access may require a Linux environment provided by the course infrastructure or campus VPN.

