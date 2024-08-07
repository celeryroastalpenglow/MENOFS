import pandas as pd

df_c_opt = pd.read_csv(snakemake.input[1],sep='\t')
c_opt = df_c_opt.level.values[0]

def change_slack(year,spatial,technology,country,optimisation,slacklevel):
	with open(snakemake.input[0], "r") as f:
		list_of_lines = f.readlines()
	list_of_lines[7] = slacklevel + "\n"
	list_of_lines[2] = str(c_opt) + "\n"
	with open("/cluster/work/projects/ec85/MENOFS_git/EU10CS/results/" + year + "_" + spatial + "_" + country + "_" + technology + "_" + optimisation + "_" + slacklevel + "/mga_parameters.dd", "w") as f:
		f.writelines(list_of_lines)

print(snakemake.wildcards.technology,snakemake.wildcards.country, snakemake.wildcards.optimisation)

change_slack(snakemake.wildcards.year,snakemake.wildcards.spatial,snakemake.wildcards.technology,snakemake.wildcards.country,snakemake.wildcards.optimisation,snakemake.wildcards.slacklevel)

