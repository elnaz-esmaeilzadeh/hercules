# System definition
# -----------------

SYSTEM = MACBOOK

# Main compilation setup
# ----------------------

#CFLAGS += -fastsse -O3
CFLAGS += -g -O0 -DDEBUG

# Compilation flags
# -----------------

# CVM and IO related:

IO_CFLAGS = -DUSECVMDB -DSCEC -DPROCPERNODE=1 -DNO_OUTPUT
IO_CPPFLAGS = -DUSECVMDB -DSCEC -DPROCPERNODE=1 -DNO_OUTPUT

# Earthquake related

SOLVE_CFLAGS = -DHALFSPACE -DBOUNDARY
SOLVE_CPPFLAGS = -DHALFSPACE -DBOUNDARY
