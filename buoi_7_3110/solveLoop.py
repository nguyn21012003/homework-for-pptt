from numpy import typing as npt


def solveLoop(AMatrix, BMatrix, dim: int, F: npt.NDArray, N: int, initInput) -> list:
    x = initInput
    listX = [x]  ### Lưu giá trị vào mảng để kiểm soát
    for i in range(1, N):
        x = F(AMatrix, BMatrix, dim, x)
        listX.append(x)
        print(x, F.__qualname__, f"tại k = {i}", "\n")

        if abs(max(listX[i]) - max(listX[i - 1])) / max(listX[i]) <= 1e-15:

            break

    return listX, i
