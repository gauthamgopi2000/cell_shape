function dx = pbc_separation(d, box_length)
dx = d - box_length * round(d / box_length);

end