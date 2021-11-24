import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("log.csv", names=[
                 "time", "primers", "mean distance", "mean number of tangent primers"])
df["time"] = df["time"] - df["time"][0]
print(df)

sns.lineplot(x="time", y="primers", data=df)
plt.show()
