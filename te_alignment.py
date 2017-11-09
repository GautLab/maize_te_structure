#The overlap function that tests of TE and RNA overlap
#this code is copied into the other file (find_te_aligns)
from random import randint

def lines_up(x1, y1, x2, y2):
    """Given a start position and an end position of 2 sequences,
    determines their alignment
    0 - No alignment
    1 - Partial alignment (even a one base alignment counts)
    2 - First totally encapsulated by second
    3 - Second totally encapsulated by first
    4 - First and second exactly match
    5 - error
    """
    if x1 == x2 and y1 == y2:
        return 4
    if x1 <= x2 and y1 >= y2:
        return 3
    if x1 >= x2 and y1 <= y2:
        return 2
    if (x2 <= x1 and x1 <= y2) or (x2 <= y1 and y1 <= y2):
        return 1
    if y1 < x2 or y2 < x1:
        return 0
    return 5

def test_lines_up(n):
    results = ["0 - No alignment",
    "1 - Partial alignment (even a one base alignment counts)",
    "2 - RNA totally encapsulated by TE",
    "3 - TE totally encapsulated by RNA",
    "4 - TE and RNA exactly match",
    "5 - error"]
    for i in range(1,n+1):
        print("TEST",i)
        print("==========")
        while (True):
            x1,y1,x2,y2=randint(1,30),randint(1,30),randint(1,30),randint(1,30)
            x1,y1=min(x1,y1),max(x1,y1)
            x2,y2=min(x2,y2),max(x2,y2)
            if x1 != y1 and x2 != y2:
                break
        print("RNA {0:2d}  {1:2d}".format(x1,y1)+" "+" "*(x1-1)+"-"*(y1-x1+1))
        print("TE  {0:2d}  {1:2d}".format(x2,y2)+" "+" "*(x2-1)+"-"*(y2-x2+1))
        print("==========")
        print(results[lines_up(x1,y1,x2,y2)])
        print()
        return
