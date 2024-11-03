{-------------------------------------------------------------------------------
- Class: CS 581
- Assignment: 3 - Abstract Syntax 
- Name: Alexis Doyle & Dhir Katre
- Date: 30 oct  2024
-------------------------------------------------------------------------------}

-- Exercise 1 Imperitive Language 

-- consider this absteact syntax of a very simple imperative languge 

type Name = String 
data Fun = Succ | Add Name 
data Stmt = Assign Name Int | Apply Fun Name | Twice Stmt 
type Prog = [Stmt]
type State = [(Name, Int)]


-- helper functions 
-- remove val 
removeVar :: Name -> State -> State
removeVar _ [] = []
removeVar x ((name, val):rest)
    | x == name = rest
    | otherwise = (name, val) : removeVar x rest


-- retrived value of a state returning 0 if not found 
getVal :: Name -> State -> Int
getVal _ [] = 0
getVal x ((name, val):rest)
    | x == name = val
    | otherwise = getVal x rest


-- define two funtions semStmt for the semantics of individul statements and semProg 
-- for the semantics of programs. semProg should call semStmt

semStmt ::  Stmt -> State -> State
semStmt (Assign x i) state = (x, i) : removeVar x state
semStmt (Apply f x) state = 
    case f of
        Succ -> (x, getVal x state + 1) : removeVar x state
        Add y -> (x, getVal x state + getVal y state) : removeVar x state
semStmt (Twice s) state = semStmt s (semStmt s state)


-- CHANGE 
semProg :: Prog -> State
semProg = foldl (flip semStmt) []

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Excersie 2 Mini Logo 

data Cmd = Pen Mode | MoveTo Int Int | Sequ Cmd Cmd
data Mode = Up | Down

type State2 = (Mode,Int,Int)                                                    -- the drawing state

type Line = (Int,Int,Int,Int)                                                   -- the points for the lines 
type Lines = [Line]

-- Semantics function for each command, maintaining state updates internally
sem :: Cmd -> State2 -> Lines
sem (Pen mode) (currentMode, x, y) = []                                         -- Only updates the pen mode, no line drawn
sem (MoveTo x y) (Up, _, _) = []                                                -- Move pen with mode Up, no line drawn
sem (MoveTo x y) (Down, x0, y0) = [(x0, y0, x, y)]                              -- Move pen with mode Down, draw line
sem (Sequ cmd1 cmd2) state =
    let lines1 = sem cmd1 state                                                 -- Process cmd1, producing its lines
        newState = updateState cmd1 state                                       -- Update state after cmd1
        lines2 = sem cmd2 newState                                              -- Process cmd2 with the new state
    in lines1 ++ lines2                                                         -- Combine lines from both commands

-- update state after each command
updateState :: Cmd -> State2 -> State2
updateState (Pen mode) (_, x, y) = (mode, x, y)
updateState (MoveTo x y) (mode, _, _) = (mode, x, y)
updateState (Sequ _ cmd2) state = updateState cmd2 state

run :: Cmd -> Lines
run cmd = sem cmd (Up, 0, 0)                                                    -- Start with pen up at (0,0)

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Excersie 3 Stack Language 

-- Consider the language S of stack programs, defined by the following grammar
-- S ::= Op | Op;S
-- Op ::= LD Int | ADD | SWAP | DUP

type S = [OP]
data OP = LD Int | ADD | SWAP | DUP

type Stack = [Int]

-- define the semantics for the stack language as SemS and define the approperate type for D
type D = Stack -> Maybe Stack

semS :: S -> D
semS [] x = Just x                                                                  -- handles empty list 
semS x s = runOP x (Just s)                                                         -- handles non-empty list

runOP :: S -> Maybe Stack -> Maybe Stack
runOP x Nothing = Nothing                                                           -- if stack nothing then error
runOP [] (Just x) = Just x                                                          -- empty list return 
runOP (x:xs) (Just s) = runOP xs (semOP x s)                                        -- recursive call to run OP

-- define an Auxilary funtion semOp 
semOP:: OP -> D
semOP (LD x) s = Just (x:s)                                                       -- push onto stack 
semOP (ADD) (x:y:xs) = Just ((x + y):xs)                                          -- add top two elements
semOP ADD _ = Nothing                                                             -- error handle for less then 2 elements 
semOP (SWAP) (x:y:xs) = Just (y:x:xs)                                             -- swap top two elements
semOP SWAP _ = Nothing                                                            -- error handle for less then 2 elements 
semOP (DUP) (x:xs) = Just (x:x:xs)                                                -- copies top element 
semOP (DUP) [] = Nothing                                                          -- fail on empty stack 
