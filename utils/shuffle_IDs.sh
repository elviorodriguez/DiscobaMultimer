#!/bin/bash

# This script shuffles the IDs in the pair list to balance the load on the GPUs.

IDs_list=$1

shuf $IDs_list > $IDs_list.shuf
