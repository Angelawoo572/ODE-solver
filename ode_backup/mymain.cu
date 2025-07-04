int main(int argc, char **argv)
{
    int nspin = 32 ;
    nBytes = 3*nspin*sizeof(float) ; 

    for(i=0,i<nspin,i++)
    {
	    ix=3*i;
	    iy=ix+1;
	    iz=iy+1;

	    if(i==0)
	    {
		    yin[ix]=0.0;
		    yin[iy]=0.0;
		    yin[iz]=1.0;
	    }
	    elseif(i == nspin)
	    {
		    yin[ix]=0.0;
		    yin[iy]=0.0;
		    yin[iz]=-1.0;
	    }
	    elseif(i < nspin/2)
	    {
		    yin[ix]=0.0;
		    yin[iy]=0.0;
		    yin[iz]=1.0;
	    }
	    elseif
	    {
		    yin[ix]=0.0;
		    yin[iy]=0.0;
		    yin[iz]=-1.0;
	    }
    }
    cudaMemcpy(y, yin, nBytes, cudaMemcpyHostToDevice);
}
