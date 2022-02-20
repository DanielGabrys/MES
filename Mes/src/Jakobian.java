
public class Jakobian
{
    double detJ;
    double[][] matrixJ=new double[2][2];
    double[][] inv_m=new double[2][2];

    Jakobian()
    {
       for(int i=0;i<2;i++)
       {
           for(int j=0;j<2;j++)
           {
            matrixJ[i][j]=0;
            inv_m[i][j]=0;
           }

       }
    }
}
