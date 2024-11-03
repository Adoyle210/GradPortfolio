{-------------------------------------------------------------------------------
- Class: CS 581
- Assignment: 3 - Abstract Syntax 
- Name: Alexis Doyle & Dhir Katre
- Date: 21 oct  2024
-------------------------------------------------------------------------------}
import Data.List (nub)   -- removes duplicate elements from a list 

-- Excerise 1 Mini Logo:

-- (a) Define the abstract sytanx for Mini Logo as a Haskell data type Cmd
-- Pen mode: can be 'up' or 'down'
data Mode = Up | Down deriving (Show)

-- Position can be either a number (Int) or a variable (String)
data Pos = Num Int | Name String deriving (Show)


data Cmd = Pen Mode                                             -- pen mode change
         | MoveTo Pos Pos                                       -- move to a position (x, y)
         | Def String [String] Cmd                              -- function definition with parameters
         | Call String [Int]                                    -- function call with arguments
         | Seq Cmd Cmd                                          -- sequence of commands
         deriving (Show)


-- (b) Write a Mini Logo function vector that draws a line from a given position (x1,y1) to a given position (x2,y2),
-- and represent the function in abstract syntax, that is, as a Haskell value of type Cmd as defined in part (a)

-- Concrete syntax (optional)
-- def vector x1 y1 x2 y3 
--  Pen Down 
--  MoveTo (x1, y1)
--  MoveTo (x2, y2)
-- Pen Up

--Vector :: Cmd 

-- Abstract syntax 
vector = Def "vector" ["x1", "y1", "x2", "y2"]
    (Seq (Pen Down)
    (Seq (MoveTo (Name "x1") (Name "y1"))
    (Seq (MoveTo (Name "x2") (Name "y2"))
        (Pen Up))))

-- OR
-- Concrete syntax (optional)
--  def vector x1 y1 x2 y3 
--  MoveTo (x1, y1)
--  Pen Down 
--  MoveTo (x2, y2)
-- Pen Up

-- -- Abstrct Syntax 
-- vector = Def "vector" ["x1", "y1", "x2", "y2"]
--     (Seq (Pen Up)
--     (Seq (MoveTo (Name "x1") (Name "y1"))
--     (Seq (Pen Down)
--     (MoveTo (Name "x2") (Name "y2")))))

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Excercise 2 Grammar Grammar 

-- Consider the following: 
-- Grammar  ::= prod; ... ; prod  
-- prod     ::= nt ::= rhs | ... | rhs 
-- rhs      ::= symbol*  (this means zero or more)
-- symbol   ::= nt | term

-- (a) Give Haskell (data) type definitions for types Grammar, Prod (production groups), RHS (right-hand side), and Symbol to represent the abstract syntax for
-- the above language. As Part of your def, use the types NT (nonterminal) and Term defined as follows

type NT = String
type Term = String

data Grammar = Grammar [Prod] deriving (Show)

data Prod = Prod NT [RHS] deriving (Show)

data RHS = RHS [Symbol] deriving (Show)

data Symbol = N NT | T Term deriving (Show)

-- (b) Consider the following for a small imperative languge Imp 
--                  cond     ::= T | not cond | (cond)
--                  stnt     ::= skip | while cond do {stnt} | stnt; stnt 
-- Represent this Grammar by a value of type Grammar defined in part (a).

imp :: Grammar 
imp = Grammar [condProd, stmtProd]

-- Auxilliary definitions 
-- Cond production 
condProd :: Prod
condProd = Prod "cond" [condRHS1, condRHS2, condRHS3]

-- T
condRHS1 :: RHS 
condRHS1 = RHS [T "T"]

-- not cond
condRHS2 :: RHS 
condRHS2 = RHS [T "not", N "cond"]

-- (cond)
condRHS3 :: RHS 
condRHS3 = RHS [T "(", N "cond", T ")"]

-- Stnt production 
stmtProd :: Prod 
stmtProd = Prod "stmt" [stmtRHS1, stmtRHS2, stmtRHS3]

stmtRHS1 :: RHS
stmtRHS1 = RHS [T "skip"]

stmtRHS2 :: RHS 
stmtRHS2 = RHS [T "while", N "cond", T "do", T "{", N "stmt", T "}"]

stmtRHS3 :: RHS
stmtRHS3 = RHS [N "stmt", T ";", N "stmt"]


-- (c) Define the following two funtions for extracting all defined nonterminals and all used terminals form a Grammar 

nonterminals :: Grammar -> [NT]
nonterminals (Grammar prods) = nub (getNTs prods)                   -- removes duplicate elements from a list 
  where
    getNTs :: [Prod] -> [NT]
    getNTs [] = []
    getNTs (Prod nt _ : rest) = nt : getNTs rest

terminals :: Grammar -> [Term]
terminals (Grammar prods) = nub (getTerms prods)                   -- removes duplicate elements from a list 
    where
        getTerms :: [Prod] -> [Term]
        getTerms [] = []
        getTerms (Prod _ rhss : rest) = getRHSTerms rhss ++ getTerms rest

        getRHSTerms :: [RHS] -> [Term]
        getRHSTerms [] = []
        getRHSTerms (RHS symbols : rest) = getSymbolTerms symbols ++ getRHSTerms rest

        getSymbolTerms :: [Symbol] -> [Term]
        getSymbolTerms [] = []
        getSymbolTerms (T t : rest) = t : getSymbolTerms rest
        getSymbolTerms (_ : rest) = getSymbolTerms rest




