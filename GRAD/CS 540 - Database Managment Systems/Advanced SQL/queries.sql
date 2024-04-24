-- Return names of every employee who works in the "Hardware", "Software", and "Research" departments. (1.5 points)
SELECT DISTINCT e.ename
FROM Emp e
JOIN Works w1 ON e.eid = w1.eid
JOIN Dept d1 ON w1.did = d1.did AND d1.dname = 'Hardware'
JOIN Works w2 ON e.eid = w2.eid
JOIN Dept d2 ON w2.did = d2.did AND d2.dname = 'Software'
JOIN Works w3 ON e.eid = w3.eid
JOIN Dept d3 ON w3.did = d3.did AND d3.dname = 'Research';


-- Return the names of every department without any employee. (1.5 points)
SELECT d.dname
FROM Dept d
LEFT JOIN Works w ON d.did = w.did
WHERE w.eid IS NULL;


-- Print the managerid of managers who manage only departments with budgets greater than $1.5 million. (1.5 points)
SELECT d.managerid
FROM Dept d
GROUP BY d.managerid
HAVING MIN(d.budge) > 1.5e6;


-- Print the name of employees whose salary is less than or equal to the salary of every employee. (1.5 points)
SELECT e.ename
FROM Emp e
WHERE e.salary <= (SELECT MIN(salary) FROM Emp);
