int mapj(int j, int nt)
{
    int i;
    if (j>=nt) {i=j-nt;}
    else if (j<0) {i=nt+j;}
    else {i=j;}
    return i;
}

