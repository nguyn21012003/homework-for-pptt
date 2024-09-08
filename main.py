import numpy as np


# Câu 1:
def converter():
    degree = float(input("Input degree in C/F: "))
    choose = input("Convert F to C or C to F: ")
    match choose:
        case "F to C":
            F = (degree - 32) * 5 / 9
            return f"{float(F):.2f}°C"
        case "C to F":
            C = degree * 9 / 5 + 32
            return f"{float(C):.2f}°F"


##########################################################


# Câu 2:
def string_to_float():
    str_num = str(input("Input a string of number: "))
    str_num = float(str_num) + 56
    return str(str_num)


##########################################################


# Câu 3:
def format_float(argument):
    formated_argument = f"{argument:.4f}"
    sqare = argument**2
    cube = argument**3
    return f"format number is {formated_argument} annd its square is {sqare:.4f} and its cube is {cube:.4f} "


##########################################################


# Câu 4:
def create_table():
    step = 2.2
    start = -10
    first_col = []
    second_col = []
    third_col = []
    for i in range(10):
        first_col.append(f"{float(start):.1f}")
        start = start + step

    for i in range(50, 70):
        if i % 2 == 0:
            second_col.append(i)

    for i in range(0, 10):
        num = i**4
        third_col.append(num)

    for i in range(len(first_col)):
        print(f" {first_col[i]: <5} | {second_col[i]} | {third_col[i]} ")


##########################################################


# Caua 5:
def find_the_smallest_num(a, b, c):
    smallest = a
    if smallest > b:
        smallest = b
    if smallest > c:
        smallest = c
    return smallest


##########################################################


# Câu 6:
def square_of_numbers():
    i = 0
    while i < 70:
        print(i**2, f"tương ứng với {i}")
        i += 1


##########################################################


# Câu 7:
def loop_range():
    for i in range(7, 14):
        number = (i + 1) ** 2
        return number


##########################################################


# câu 8
def complete_code():
    i = 0
    while True:
        for i in range(10, 80):
            print((i + 1) ** 2, "tương ứng với", i + 1)
            i += 1
            if i == 80:
                return False


##########################################################


# câu 9:
def create_list_number(n):
    n = int(n)
    list_nums = []
    for i in range(n):
        list_nums.append(i)

    with open("list_number.txt", "w", newline="") as writefile:
        for row in list_nums:
            pass


def main():
    """
    Các bài tập được chia thành từng block code, kiểm tra từng câu bằng cách bỏ kí tự hashtag "#" ở câu đó ở bên dưới
    """
    # Câu 1:
    # print(converter())
    ##########################################################
    # Câu 2:
    # print(string_to_float())
    ##########################################################
    # Câu 3:
    # print(format_float(27847.91284582))
    ##########################################################
    # Câu 4:
    # create_table()
    ##########################################################
    # Câu 5:
    # smallest = find_the_smallest_num(1, 2, 0.1)
    # print(smallest)
    ##########################################################
    # Câu 6
    # square_of_numbers()
    ##########################################################
    # câu 7
    # print(loop_range())
    ###########################################################
    # câu 8
    complete_code()
    ##########################################################


if __name__ == "__main__":
    main()
