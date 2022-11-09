from findallfiles import findAllfiles
from os import remove
pa = 'D:\\Programming\\User'
findAllfiles(pa)

f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'r')

# print(li)
try:
    for i in f.readlines():
        suffix = i.split(sep='.')[-1][:-1]
        # print(i, i[-3:])
        # print('|{}|'.format(suffix), len(suffix))
        if suffix in ['spec']:
            print(i)
            remove(i[:-1])
except PermissionError:
    pass
f.close()
# print(findAllfiles())
