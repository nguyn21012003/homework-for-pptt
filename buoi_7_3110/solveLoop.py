from numpy import typing as npt
from tqdm import tqdm


def solveLoop(AMatrix: npt.NDArray, BMatrix: npt.NDArray, dim: int, F: npt.NDArray, N: int, initInput) -> list:

    eigenFunction = initInput
    listEigenFunction = [eigenFunction]  ### Lưu giá trị vào mảng để kiểm soát
    for i in tqdm(range(1, N), desc=f"{F.__qualname__}", unit="step"):

        eigenFunction = F(AMatrix, BMatrix, dim, eigenFunction)
        listEigenFunction.append(eigenFunction)
        # print(eigenFunction, F.__qualname__, f"tại k = {i}", "\n")

    # s = np.array(listEigenFunction)

    return eigenFunction, i
