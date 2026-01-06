import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the datasets
df_pooled = pd.read_csv('optimal_timeslots_pooled.csv')
df_separate = pd.read_csv('optimal_timeslots_separate.csv')

# Set up the plotting style
sns.set(style="whitegrid")

# Create a figure with subplots
fig = plt.figure(figsize=(16, 12))

# 1. Boxplot comparison of Waiting Days
ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
# Combine for plotting
df_combined = pd.concat([df_pooled, df_separate])
sns.boxplot(x='policy', y='waiting_days', data=df_combined, ax=ax1, palette="Set2")
ax1.set_title('Distribution of Mean Waiting Days: Separate against Pooled Policy', fontsize=16)
ax1.set_ylabel('Mean Waiting Days', fontsize=12)
ax1.set_xlabel('Policy', fontsize=12)

# 2. Heatmap for Separate Policy
ax2 = plt.subplot2grid((2, 2), (1, 0))
pivot_sep = df_separate.pivot(index="ts2", columns="ts1", values="waiting_days")
sns.heatmap(pivot_sep, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax2, cbar_kws={'label': 'Waiting Days'})
ax2.set_title('Separate Policy: Waiting Days by Timeslot Size', fontsize=14)
ax2.set_xlabel('Type 1 Timeslot (min)', fontsize=10)
ax2.set_ylabel('Type 2 Timeslot (min)', fontsize=10)

# 3. Heatmap for Pooled Policy
ax3 = plt.subplot2grid((2, 2), (1, 1))
pivot_pool = df_pooled.pivot(index="ts2", columns="ts1", values="waiting_days")
sns.heatmap(pivot_pool, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax3, cbar_kws={'label': 'Waiting Days'})
ax3.set_title('Pooled Policy: Waiting Days by Timeslot Size', fontsize=14)
ax3.set_xlabel('Type 1 Timeslot (min)', fontsize=10)
ax3.set_ylabel('Type 2 Timeslot (min)', fontsize=10)

plt.tight_layout()
plt.savefig('policy_comparison_waiting_times.png')