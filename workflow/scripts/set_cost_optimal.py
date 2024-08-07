import sqlite3
import pandas as pd

path = "/cluster/work/projects/ec85/MENOFS_git/models/2010/region"
db_file = path + '/results.db'
costs=[]
con = sqlite3.connect(db_file)
(
    costs
    .append(
        (
            round(pd
            .read_sql_query("SELECT * from scalarvariables", con),3)
        )
    )
)
con.close()
df_costs = pd.concat(costs)

df_costs.to_csv(snakemake.output.cost_opt,sep='\t')