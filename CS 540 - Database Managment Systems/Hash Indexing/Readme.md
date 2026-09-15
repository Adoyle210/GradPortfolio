# Hash Indexing

A static hash index built on top of a page-organized storage layer, giving near-constant-time employee lookups instead of scanning pages sequentially.

## Summary
This extends the page/slot-directory storage design from the **External Storage Management** project by adding a hash index: each `Employee` record's ID is hashed to a bucket, and a page directory maps each bucket to the page on disk holding its records. Looking up an employee by ID becomes "hash the ID, jump straight to its page" instead of scanning every page in the file.

## Tech Stack
- C++
- Binary file I/O
- Static hashing with a page directory and overflow-page support

## How it works
- **Records and pages** reuse the same slotted-page format as the storage management project — 4KB pages, variable-length records, per-page slot directories.
- **Hashing**: each record's ID is run through `compute_hash_value()` to determine its bucket.
- **Page directory**: an in-memory table mapping each hash bucket to the offset of its page in the index file (`EmployeeIndex`). The first time a bucket is used, a new page is allocated for it; later records with the same bucket are appended to that existing page.
- **Overflow pages**: each page reserves a pointer for chaining to an overflow page if a bucket's records outgrow a single 4KB page.
- **Lookup**: hashes the requested ID, looks up its page directly via the page directory, and scans only that one page (plus any overflow pages) for a match — no full-file scan required.

## Build & Run
```sh
g++ -std=c++11 main.cpp -o main.out
./main.out
```
On startup, the program reads `Employee.csv`, builds the hash index in `EmployeeIndex`, then prompts for employee IDs to look up (type `exit` to quit).

## Notes
- This project is best presented alongside **External Storage Management** — together they show the same on-disk page format used for both a straightforward sequential store and an indexed one, and the concrete lookup-speed benefit indexing adds.
