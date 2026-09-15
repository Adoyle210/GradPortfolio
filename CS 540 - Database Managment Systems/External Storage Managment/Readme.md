# External Storage Management

A C++ program that implements page-organized storage for variable-length records on disk, following the classic slotted-page layout used by real database systems.

## Summary
Given an `Employee(id, name, bio, manager_id)` relation as a CSV file, this program serializes the records into fixed-size 4KB pages — each with a slot directory tracking record offsets and sizes — and writes those pages to a binary data file. It then supports looking up an employee by ID directly from disk, reading only a few pages into memory at a time rather than loading the whole file.

## Tech Stack
- C++
- Binary file I/O
- Slotted-page storage design (per Ramakrishnan & Gehrke's *Database Management Systems*, Ch. 9)

## How it works
- **Records** are variable-length: `name` and `bio` are stored as length-prefixed strings rather than fixed-size buffers, so each record only takes as much space as it needs.
- **Pages** are exactly 4096 bytes. Each page has a slot directory — a list of `(offset, size)` pairs — that records where each record lives inside the page, so a record can be located without scanning the whole page byte-by-byte.
- **Storage manager** keeps at most 3 pages in memory at once while building the data file, flushing full pages to disk (`EmployeeRelation.dat`) as it goes, and reads pages back the same way — 3 at a time — while searching.
- **Lookup** takes an employee ID from the console, walks the on-disk pages using their slot directories, and prints the matching record (repeatable — you can search multiple IDs in one run without restarting).

## Build & Run
```sh
g++ -std=c++11 main.cpp -o main.out
./main.out
```
On startup, the program reads `Employee.csv` from the current directory and builds `EmployeeRelation.dat`. It then prompts for employee IDs to search for.

## Notes
- Building on this design, see the companion **Hash Indexing** project, which adds a hash index directly on top of this same page format for faster lookups.
