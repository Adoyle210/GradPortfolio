-- Return names of every employee who works in the "Hardware", "Software", and "Research" departments. (1.5 points)
SELECT E.ENAME FROM EMP E, DEPT D1, DEPT D2, WORKS W1, WORKS W2
WHERE
E.EID = W1.EID AND W1.DID = D1.DID AND D1.DNAME ='HARDWARE' AND
E.EID = W2.EID AND W2.DID = D2.DID AND D2.DNAME ='SOFTWARE' AND
E.EID = W2.EID AND W2.DID = D2.DID AND D2.DNAME ='RESEARCH'
OR
SELECT E.ENAME FROM EMP E, DEPT D, WORKS W
WHERE
E.EID = W.EID AND W.DID = D.DID AND D.DNAME ='HARDWARE'
INTERSECT
SELECT E1.ENAME FROM EMP E1, DEPT D1, WORKS W1
WHERE
E1.EID = W1.EID AND W1.DID = D1.DID AND D1.DNAME ='SOFTWARE'  --- fix to add reasearch 


-- Return the names of every department without any employee. (1.5 points)


-- Print the managerid of managers who manage only departments with budgets greater than $1.5 million. (1.5 points)
SELECT D.MANAGERID FROM DEPT D WHERE D.BUDGET > 1500000
EXCEPT 
SELECT D1.MANAGERID FROM DEPT D1 WHERE D1.BUDGET <= 1500000

-- Print the name of employees whose salary is less than or equal to the salary of every employee. (1.5 points)