# -------------------------------------------------------------------------------
# Class: CS 581
# Assignment: 1 - warm up
# Name: Alexis Doyle 
# Date: 25 Sep 2024
# -------------------------------------------------------------------------------
# (A) Languge: Python 
# (B) Directions: 
#     Please just run the python file main, I have already hard coded in inputs that match the 
#     samples for ease of grading. Each exercise is seprated in its own def below with comments. 
#     Everything is run in main if you want to edit any of the inputs. Thanks :) 



# Exercise 1: Counting Specific Numbers
# count the numbers of fives or threes in the list. 
# -------------------------------------------------------------------------------
def count(x):
    print(f'input: {x}')
    num = 0

    for i in x:
        if i == 5 or i == 3:
            num += 1

    print (f'output: {num}')
    print (' ')

# Exercise 2: Longest Ascending Subsequence 
# Determine the longest subsequence of consecutive non-decreasing intergers. 
# In the case of multiple solutions, any longest subsequence is acceptable. 
# -------------------------------------------------------------------------------
def longest (x):
    print(f'input: {x}')
    nums = []  # longest asending subsequence 
    temp = []  # temp list to find  LAS

    nlen = 0   # length of nums
    tlen = 0   # length of temp 
    prev = None   #  the previous number in the list 

    for i in x :
        if prev is None or i >= prev:                 
            tlen += 1 
            temp.append(i)
            prev = i
        else:
            if nlen < tlen:
                # store current longest
                nlen = tlen 
                nums = temp.copy()
            # start new temp list  
            temp.clear()
            temp = [i]
            tlen = 1
            prev = i

    #final check 
    if nlen < tlen:
        nums = temp.copy()

    print (f'output: {nums}')
    print (' ')




# Exercise 3: Reverse Polish Notation
# In Reverse polish notatior (RPN) expressions are written as a sequence of numbers and operations, 
# without the need for parathesis. In particullar, operations follow thier operands. 
# -------------------------------------------------------------------------------
def RPN (x):
    print(f'input: {x}')                                #["2","3","+","4","*"] => ["20"]           2 3 + 
    stack = []

    for i in x :

        #case for adding opperation
        if i == '+':
           #popping values off stack 
           b = stack.pop()
           a = stack.pop()

           temp = a + b

           #adding sol to stack 
           stack.append(temp)

        #case for multiplying opperation
        elif i == '*' :
           #popping values off stack 
           b = stack.pop()
           a = stack.pop()

           temp = a * b

           #adding sol to stack 
           stack.append(temp)

        #case for the square operation
        elif i == 'sqr' :
           #popping values off stack 
           b = stack.pop()

           temp = b * b 

           #adding sol to stack 
           stack.append(temp)
        
        else:
            #adding numbers 
            stack.append(int(i)) 

        print(f'   stack: {stack}')

    print (f'output: {stack}')
    print (' ')
  


# main 
# used to run the program
# -------------------------------------------------------------------------------
def main ():
    print('Welcome to my Homework warm up! \n  Running exercise 1 ... \n ')

    # Exercise 1:
    count([2,3,5,1,3,3,7])
    count([2,7,8,4])
    count([3,8,2,4,3])

    print('Exercise 1 Completed :) \n  Running exercise 2 ... \n')

    # Exercise 2:
    longest([2,3,5,1,3,3,7])
    longest([9,7,8,4])
    longest([8,4,3,1])

    print('Exercise 2 Completed :) \n  Running exercise 3 ... \n')

    # Exercise 3:
    RPN(["2", "3", "+", "4", "*"])
    RPN(["2", "3", "4", "+", "*"])
    RPN(["2", "3", "2", "sqr","+", "*"])
    RPN(["5", "1", "2", "+","4", "*","+", "1", "+"])  #not an example just a test

    print('All the exercises are now complete :) ')



main()







