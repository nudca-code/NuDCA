import sys
from pathlib import Path
sys.path.insert(1, str(Path.cwd().resolve().parent))

import matplotlib.pyplot as plt
from nudca.decay_database import load_decay_database

decay_database = load_decay_database(data_source='ENDF-B-VIII.1_decay')

fig = decay_database.plot_nuclear_chart()

plt.tight_layout()
plt.show()