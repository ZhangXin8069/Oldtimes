import pandas as pd

#
# data = {
#     "calories": [420, 380, 390],
#     "duration": [50, 40, 45]
# }
#
# df = pd.DataFrame(data, index=["day1", "day2", "day3"])

# 指定索引
# print(df.loc["day2"])
df = pd.read_excel(r"D:\Programming\User\txt\12月账单(2)(1)(1).xlsx")
print(df.to_string())
df.to_csv('zx.txt', index=False, encoding='gbk')
