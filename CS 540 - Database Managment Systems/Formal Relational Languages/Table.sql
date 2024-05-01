-- set up all the tables in the database 

DROP TABLE IF EXISTS Emp;
 
CREATE TABLE Emp (            
    eid int NOT NULL AUTO_INCREMENT,    
    estring varchar(25) NOT NULL,
    age int,
    salary real,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS Works;
 
CREATE TABLE Works (            
    eid int NOT NULL AUTO_INCREMENT,    
    did int,
    pc_time int,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS Dept;
 
CREATE TABLE Dept (            
    did int NOT NULL AUTO_INCREMENT,    
    dstring varchar(25) NOT NULL,
    budget real,
    managerid int,
    PRIMARY KEY (id)
);
