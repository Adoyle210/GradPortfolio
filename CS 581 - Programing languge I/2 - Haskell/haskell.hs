{-------------------------------------------------------------------------------
- Class: CS 581
- Assignment: 2 - haskell
- Name: Alexis Doyle & Dhir Katre
- Date: 8 oct  2024
-------------------------------------------------------------------------------}

-- Excerise 1 - Programming with lists 

type Bag = [(Int,Int)]

list = [2,3,3,5,7,7,7,8]

bag1 :: Bag
bag1 = [(5,1),(7,3),(2,1),(3,2),(8,1)]

-- (a) Define the function ins that inserts an element into a multiset 
ins :: Int -> Bag -> Bag
ins x [] = [(x, 1)]                                                                  --  if bag empty make new bag         
ins x ((y, n):ys)                                                                    --  match non-empty bags
    | x == y       = (y, n + 1) : ys                                                 -- if match increment 
    | otherwise    = (y, n) : ins x ys                                               -- if no match add 

-- (b) Define the function del that removes a single element from a multiset
del :: Int -> Bag -> Bag
del a []     = []                                                                   -- del from empty is empty 
del a ((x,n):xs)                                                                    -- decrease count of x and constuct new list removing items 0 or less
            | a == x     = if n > 1 then (x, n-1):xs else xs          
            | otherwise  = (x,n) : del a xs                                         -- mo match keep checking 

-- (c) Define a function subbag that determines whether its first argument bag is contained in the second.
subbag :: Bag -> Bag -> Bool 
subbag [] _ = True                                                                  -- if first bag empty then subbag 
subbag _ [] = False                                                                 -- if secound bag is empty then false 
subbag (x:xs) ys = check x ys && subbag xs ys                                       -- recursive check 


check :: (Int, Int) -> Bag -> Bool
check (x, n) []  = False                                                            -- second bag is empty and no match false
check (x, n) ((y, m):ys)                                                            -- checking first bag againt second 
        | x == y   = n <= m                                                         -- found match :)
        | otherwise = check (x, n) ys                                               -- keep looking 

-- (d) Define a function isSet that tests whether a bag is actually a set, which is the case when each element occurs only once.
isSet :: Bag -> Bool 
isSet [] = True                                                                     -- empty is true
isSet ((_, count): xs) = count == 1 && isSet xs                                     -- check each time is once 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Exercise 2 - Counting Specific Numbers (only 5s and 3s)
count :: [Int] -> Int 
count []  =  0                                                                      -- if empty then 0
count (x:xs)
    | x == 5 || x == 3 = 1 + count xs                                               -- only add 1 for 5 and 3 
    | otherwise            = count xs                                               -- otherwise just keep looking 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Exercise 3 - Longest Ascending Subsequence 

longestAscending :: [Int] -> [Int]
longestAscending [] = []                                                                -- empty is empty
longestAscending (x:xs) = findLongest xs [x] x []                                       -- Start with the first element

append :: [Int] -> [Int] -> [Int]
append []       ys = ys
append (x:xs)   ys = x : append xs ys

-- Recursive function to find the longest ascending subsequence
-- [Int] -> Remaining list                        
-- [Int] -> Current longest Subsequence           -> cL
-- Int   -> Last element of current Subsequence   -> lE
-- [Int] -> Longest sequence found so far         -> lF
-- [Int] -> return type

findLongest :: [Int] -> [Int] -> Int -> [Int] -> [Int]
findLongest [] cL _ lF
    | len cL > len lF = cL                                                          -- returns current if it's longer
    | otherwise = lF                                                                -- returns the longest Found

findLongest (y:ys) cL lE lF 
    | y >= lE   = findLongest ys (append cL [y]) y lF                               -- Continue current seq.
    | otherwise = findLongest ys [y] y (if len cL > len lF then cL else lF)         -- Start new sequence

-- saving while number is smaller 
ascendingFrom :: Int -> [Int] -> [Int]
ascendingFrom x xs = x : takeWhile (> x) xs                                         -- docs on takeWhile : http://www.zvon.org/other/haskell/Outputprelude/takeWhile_f.html

-- length of a list 
len :: [Int] -> Int 
len [] = 0 
len (x:xs) = 1 + len xs 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Exercise 4 - Reverse Polish Notation 
type Stack = [Int]                                                                    -- define a stack as list of ints 

rpn :: [String] -> [String]
rpn [] = [] 
rpn input = [show (head (processList [] input))]                                   

processList :: Stack -> [String] -> Stack 
processList stack [] = stack 
processList stack (item:remaining)
    | item == "+" = processList (add stack) remaining                                  -- if + add 
    | item == "*" = processList (multiply stack) remaining                             -- if * multiply 
    | item == "sqr" = processList (square stack) remaining                             -- if sqr square 
    | otherwise =  processList ((read item) : stack) remaining                         -- if its not an operation then process as num

-- operations: 
add :: Stack -> Stack 
add (y:x:remaining) = (x + y) : remaining                                             -- pop two nums, add then push 
add stack = stack                                                                     -- ignore if not two nums, to handle ["3", "+"] -> ["3"]

multiply :: Stack -> Stack
multiply (y:x:remaining) = (x * y) : remaining                                        -- Pop two nums, multiply then push 
multiply stack = stack 

square :: Stack -> Stack
square (x:remaining) = (x * x) : remaining                                           -- square then push 
