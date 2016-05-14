import timeit
from karta.vector.table import Table

def test_creation_single():
    print(timeit.timeit(stmt="Table(3.1415)", setup="from karta.vector.table import Table"))
    return

def test_creation_list():
    print(timeit.timeit(stmt="Table([1,2,3,4,5])", setup="from karta.vector.table import Table"))
    return

def test_creation_dict():
    print(timeit.timeit(stmt="Table({'value1':[1,2,3,4,5], 'value2':[6,7,8,9,10]})",
                        setup="from karta.vector.table import Table"))
    return

def test_point_indexing():
    print(timeit.timeit(stmt="m[2]", setup="from karta.vector.table import Table; m = Table([1,2,3,4,5])"))
    return


def test_field_indexing():
    print(timeit.timeit(stmt="m.getfield('b')",
                        setup="from karta.vector.table import Table;"
                              "m = Table({'a': [1,2,3,4,5], 'b':[6,7,8,9,10]})"))
    return

if __name__ == "__main__":
    print("single value", end="\n\t")
    test_creation_single()
    print("list value", end="\n\t")
    test_creation_list()
    print("dict value", end="\n\t")
    test_creation_dict()
    print("index point from list", end="\n\t")
    test_point_indexing()
    print("index field", end="\n\t")
    test_field_indexing()

