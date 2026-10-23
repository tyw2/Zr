#!/usr/bin/env python3


def generate_summary_py():
    return '''#!/usr/bin/env python3
"""汇总 EOS 与 cohesive 已完成的结果；不修改原始任务文件。"""
import csv,json
from pathlib import Path
root=Path(__file__).resolve().parent
all_eos=[]
all_cohesive=[]
all_elastic=[]
all_formation=[]
for d in sorted(root.glob('seed_*')):
    if not d.is_dir(): continue
    for e in sorted(d.glob('EOS-*/eos_result.csv')):
        with e.open() as f:
            for row in csv.DictReader(f):
                all_eos.append(dict(seed=d.name,phase=e.parent.name[4:],**row))
    for e in sorted(d.glob('Cohesive-*/cohesive_result.json')):
        row=json.loads(e.read_text())
        all_cohesive.append(dict(seed=d.name,**row))
    for e in sorted(d.glob('Formation-*/formation_result.json')):
        data=json.loads(e.read_text())
        all_formation.append(dict(seed=d.name,**data))
    for e in sorted(d.glob('Elastic-*/elastic_result.json')):
        data=json.loads(e.read_text())
        data={k:v for k,v in data.items() if not isinstance(v, (dict,list))}
        all_elastic.append(dict(seed=d.name,**data))
for fn,rows in [('eos_summary.csv',all_eos),('cohesive_summary.csv',all_cohesive),
                ('elastic_summary.csv',all_elastic), ('formation_summary.csv',all_formation)]:
    if rows:
        with (root/fn).open('w',newline='') as f:
            fields=list(rows[0])
            for row in rows[1:]:
                for field in row:
                    if field not in fields:
                        fields.append(field)
            w=csv.DictWriter(f,fieldnames=fields)
            w.writeheader(); w.writerows(rows)
        print(f'[OK] {root/fn}: {len(rows)} rows')
    else:
        print(f'[info] No finished data for {fn}')
'''
