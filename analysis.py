import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# filenames = ["data/20L-9ED-genetic.csv", "data/25L-9ED-genetic.csv", "data/30L-9ED-genetic.csv"]
filenames = ["data/20L-8ED-genetic.csv", "data/20L-9ED-genetic.csv", "data/20L-10ED-genetic.csv"]
# filenames = ["data/20L-9ED-genetic.csv", "data/20L-9ED-one-by-one.csv"]

dfs = []
for filename in filenames:
    df = pd.read_csv(filename, names=[
                     "time", "primers", "mean distance", "mean number of tangent primers"])
    df["time"] = df["time"] - df["time"][0]
    dfs.append(df)
df = pd.concat(dfs,
                keys=filenames,
                names=["file", "data point"],
                join="inner")
df = df.reset_index()
print(df)

sns.lineplot(x="time", y="primers", data=df, hue="file")
plt.show()
