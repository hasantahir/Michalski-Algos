function indx = indexx(n,arrain)
% !***********************************************************************
% !
% !     indexx - takes an array and produces a set of indices such
% !              that arrin(indx(j)) is in ascending order j=1,2,..,n
% !
% !     input data:
% !              n number of elements in the array to be ordered
% !          arrin r*8 array which is to placed in ascending order.
% !
% !     output data:
% !            indx a set of indices for ascending indices
% !
% !     notes:
% !
% !     this routine is taken from the book numerical receipes by
% !     press, flannery, teukolsky and vetterling chapter 8 p. 233.
% !     isbn 0-521-30811-9 pub. cambridge university press (1986)
% !     qa297.n866
% !
% !     this routine has been adapted by charles j gillan for use
% !     in this code.
% !
% !***********************************************************************
% !
vsmall = 1e-20;
% !---- screen for bad input data
% !
if (n <= 0)
    error('error. \ninput unspecified');
end
% !
% !---- initialize the index array with consecutive integers
% !
indx = 1 : n;
% !
% !---- we'd better remove the special case of just being called
% !     with one element! there is no more work to do in that case.
% !
% !     we will have set
% !
% !        indx(1) = 1
% !
% !     above, which is all we need to do.
% !
if(n == 1)
    return;
end
% !
% !---- ok, we get here with at least two elements.
% !
l = round(n/2 + 1);
ir = n;
% !
% !---- from here on the algorithm is heapsort with indirect addressing
% !     through indx in all references to arrin
% !
% 10  continue
% !
if(l > 1)
    l = l - 1;
    indxt = indx(l);
    q = arrain(indxt);
else
    indxt = indx(ir);
    q = arrain(indxt);
    indx(ir) = indx(1);
    ir = ir - 1;
    if(ir == 1)
        indx(1) = indxt;
        return;
    end
end
% !
i = l;
j = l+l;
% !
% 20  continue
% !
if(j <= ir)
    if(j < ir)
        if( arrin(indx(j)) < arrin(indx(j+1)) + vsmall )
            j = j + 1;
        end
        if( q < arrin(indx(j)) + vsmall )
            indx(i) = indx(j);
            i = j;
            j = j+j;
        else
            j = ir + 1;
        end
    end
    indx(i) = indxt;
end
