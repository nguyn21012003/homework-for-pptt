import sys

from random import randint


# câu 1: Write a code that converts temperature in 𝐹 [Farenheit] to 𝐶 [Celsius], and vice versa.
def converter():
    degree = input("Input degree: ")
    choose = input("F to C or C to F: ")
    match choose:
        case "F to C":
            C = (float(degree) - 32) * 5 / 9
            print(f"{C:.1f}°C")
        case "C to F":
            F = float(degree) * 9 / 5 + 32
            print(f"{F:.1f}°F")


############################################################


# câu 2: Convert the string "123.4" into a float, add 56 to it, and then convert it back into a string
def conver_string(arg):
    arg = float(arg) + 56
    return str(arg)


############################################################


# câu 3: Format the float 27847.91284582, its square, and its cube, using scientific notation with 2, 3, and 4 decimal places respectively.
def format_float_square_cube(arg):
    formatted_str = f"{arg:.4f}"
    square = arg**2
    cube = arg**3
    return f"Format number is {formatted_str} and its square is {square:.4f}, its cube is {cube:.4f}"


############################################################


# câu 4: Use string formatting and sign control to create a table of three columns of numbers. The first column should include 10 numbers starting from −10 in steps of 2.2, the second column should include all odd numbers between 50 and 70, and the third column should include all powers of 4 of the integers 0−9. Use any delimiter you wish to visually separate the three columns.”
def create_table():
    step = 2.2
    start = -10
    first_col = []
    second_col = []
    third_col = []

    for i in range(10):
        first_col.append(f"{float(start):.1f}")
        start += step

    for i in range(50, 70):
        if i % 2 == 0:
            second_col.append(i)

    for i in range(10):
        third_col.append(i**4)

    for i in range(len(first_col)):
        print(f"{first_col[i]:<5} | {second_col[i]} | {third_col[i]}")


############################################################


# câu 5: Find the smallest of three real numbers a, b, and c.
def smallest_num(a, b, c):
    smallest = a
    if smallest > b:
        smallest = b
    if smallest > c:
        smallest = c

    return smallest


############################################################


# cau 6: Print out all square numbers less than 70.
def while_loop(n):
    i = 0
    square_nums = []
    while i < n:
        sq = i**2
        square_nums.append({i: sq})
        i += 1

    return square_nums


############################################################


# câu 7: Calculate the squares of all numbers between 7 and 14.
def for_range_loop(start: int, end: int) -> list:
    """Calculate the squares of all numbers between 7 and 14.


    Args:
        start (int): start index, start with 7
        end (int): end index, end with 14

    Returns:
        sqs ( list ): contain dicts in list, a dict {i: square of i_th}
    """
    sqs = []
    for i in range(int(start), int(end)):
        if 7 < i < 14:
            sq = (i) ** 2
            sqs.append({i: sq})

    return sqs


############################################################


# câu 8: Complete the following code to print all square numbers between 10 and 80.
def complete_code(start, end):
    i = 0
    sqs = []
    while True:
        for i in range(int(start), int(end)):
            sq = i**2
            sqs.append({i: sq})
            i += 1
            if i == 80:
                print(sqs)
                return False


############################################################


# câu 9: Create a list containing the squares of all integers from 0 to 100 that are divisible by 3 and 7
def list_numb(n):
    n = int(n)
    list_of_int = []
    for i in range(n):
        square = i**2
        if square % 3 == 0 and square % 7 == 0:
            list_of_int.append(square)

    return list_of_int


############################################################


# câu 10: Write a Python program that reads in a one-column list of numbers, and calculates the sum.
def calculates_the_sum():
    list_of_numbs = []
    i = 0
    randome_int = randint(1, 4)
    while i <= randome_int:
        list_of_numbs.append(i)
        i += 1
    sum_col = 0
    for col in range(len(list_of_numbs)):
        sum_col += col

    print(list_of_numbs, end="\n")
    print(sum_col)


############################################################


# câu 11: Write a Python program that calculates the 6 times table, and outputs it to a file 6times.dat
def times_table(x, y):
    if int(x) == 0:
        sys.exit("Bảng cửu chương không có số 0!")
    elif int(x) != 1:
        sys.exit("Invalid input!")
    with open("6times.dat", "w", newline="") as writefile:
        x, y = int(x), int(y)
        muls = []
        while x <= 10:
            mul = x * y
            muls.append(mul)
            x += 1

        for row in range(len(muls)):
            writefile.write(f"{row + 1} x {y} = {muls[row]}\n")


############################################################


def main():
    """
    Các câu nhỏ của bài tập 3 được viết trong 1 block code, để kiểm tra từng câu nhỏ, ta bỏ đi dấu hashtag "#" ở mỗi câu
    """
    # câu 1
    # converter()
    ############################################################
    # câu 2
    # print(conver_string("123.4"))
    ############################################################
    # câu 3
    # print(format_float_square_cube(27847.91284582))
    ############################################################
    # câu 4
    # create_table()
    ############################################################
    # câu 5
    # smallest = smallest_num(1.0001, 4.2, 10)
    # print("The number smallest is", smallest)
    ############################################################
    # câu 6
    # print(while_loop(70))
    ############################################################
    # câu 7
    # print(for_range_loop(7, 14))
    ############################################################
    # câu 8
    # complete_code(10, 80)
    ############################################################
    # câu 9
    # print(list_numb(100))
    ############################################################
    # câu 10
    # calculates_the_sum()
    # câu 11
    # times_table(1, 6)


if __name__ == "__main__":
    main()
