nodes = dict()


# Prints out the triangle nodes if they are valid.
def triangle(x1, y1, x2, y2, x3, y3):
    if check(x1, y1) and check(x2, y2) and check(x3, y3):
        print nodes[(x1, y1)], nodes[(x2, y2)], nodes[(x3, y3)], 0.0


# Prints out the boundary nodes if they are valid.
def boundary(x, y):
    if x == 0 or y == 0:
        print nodes[(x, y)], 0.0
    elif x == 3 and y >= 4:
        print nodes[(x, y)], 15.0
    elif x >= 4 and y == 4:
        print nodes[(x, y)], 15.0


# Checks if the given (x, y) pair is actually a node.
def check(x, y):
    return 0 <= x < 6 and 0 <= y < 6 and not (x >= 4 and y == 5)


if __name__ == '__main__':
    # file1.dat
    index = 1
    for x in range(6):
        for y in range(6):
            if x >= 4 and y == 5:
                continue
            nodes[(x, y)] = index
            print index, 0.02 * x, 0.02 * y
            index += 1

    # file2.dat
    for x in range(6):
        for y in range(6):
            if x >= 4 and y == 5:
                continue
            triangle(x, y, x + 1, y, x, y + 1)
            triangle(x, y, x - 1, y, x, y - 1)

    # file3.dat
    for x in range(6):
        for y in range(6):
            if x >= 4 and y == 5:
                continue
            boundary(x, y)
