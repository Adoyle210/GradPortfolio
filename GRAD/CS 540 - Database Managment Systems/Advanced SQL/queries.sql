-- Return names of every employee who works in the "Hardware", "Software", and "Research" departments. (1.5 points)
SELECT DISTINCT e.ename
FROM emp e
JOIN works w1 ON e.eid = w1.eid
JOIN dept d1 ON w1.did = d1.did AND d1.dname = 'Hardware'
JOIN works w2 ON e.eid = w2.eid
JOIN dept d2 ON w2.did = d2.did AND d2.dname = 'Software'
JOIN works w3 ON e.eid = w3.eid
JOIN dept d3 ON w3.did = d3.did AND d3.dname = 'Research';


-- Return the names of every department without any employee. (1.5 points)
SELECT d.dname
FROM dept d
LEFT JOIN works w ON d.did = w.did
WHERE w.eid IS NULL;


-- Print the managerid of managers who manage only departments with budgets greater than $1.5 million. (1.5 points)
SELECT d.managerid
FROM dept d
GROUP BY d.managerid
HAVING MIN(d.budget) > 1.5e6;


-- Print the name of employees whose salary is less than or equal to the salary of every employee. (1.5 points)
SELECT e.ename
FROM Emp e
WHERE e.salary <= (SELECT MIN(salary) FROM Emp);



-- Print the enames of managers who manage the departments with the largest budget (0.5 point).

SELECT e.ename
FROM emp e
JOIN dept d ON e.eid = d.managerid
WHERE d.budget = (
    SELECT MAX(budget) 
    FROM dept
);


-- Print the name of every department and the average salary of the employees of that department. The department must have a budget more than or equal to $50. (0.5 point)

SELECT d.dname, AVG(e.salary) AS avg_salary
FROM dept d
JOIN works w ON d.did = w.did
JOIN emp e ON w.eid = e.eid
WHERE d.budget >= 50
GROUP BY d.did, d.dname;

-- Print the managerids of managers who control the largest amount of total budget. As an example, if a manager manages two departments, the amount of total budget for him/her will be the sum of the budgets of the two departments. We want to find managers that have max total budget. (1 point)

SELECT d.managerid
FROM dept d
GROUP BY d.managerid
ORDER BY SUM(d.budget) DESC
LIMIT 1;


-- Print the name of every employee who works only in the ”Hardware” department. (1 point)

SELECT e.ename
FROM emp e
JOIN works w ON e.eid = w.eid
JOIN dept d ON w.did = d.did
WHERE d.dname = 'Hardware'
GROUP BY e.eid, e.ename
HAVING COUNT(*) = 1;


