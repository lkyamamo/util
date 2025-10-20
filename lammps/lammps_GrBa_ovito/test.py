import pandas as pd

a = pd.DataFrame({"a":[1,2,3,4,5], "b":[6,7,8,9,10], "c":[11,12,13,14,15]})

bin_size = 1/3
names = a.columns
names = [name for name in names if name != "c"]

print(names)

prefactors = (a[names]).sum(axis=0)

for name, factor in zip(names, prefactors):
    a[name] = a[name]/factor

print(a)
