import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Grid
{
    double H; //height
    double B; //length
    int nH;
    int nB;
    int nE; //number of elements
    int nN; //number of nodes

    Element[] elements; //array of elements
    Node[] nodes; //array of nodes

    double[][] H_global;
    double[][] Hbc_global;
    double[] P_global;
    double[][] C_global;


    Grid(PointSchema p)
    {
        this.H=Data.H;
        this.B=Data.B;
        this.nH=Data.nH;
        this.nB=Data.nB;
        this.nE = (nB - 1) * (nH - 1);
        this.nN = nB * nH;

        H_global = new double[nN][nN];
        Hbc_global= new double[nN][nN];
        C_global = new double[nN][nN];
        P_global = new double[nN];

        elements = new Element[nE];
        nodes = new Node[nN];

        init();

        double dx = B / ((double)nB - 1);
        double dy = H / ((double)nH - 1);

        int row = 0;
        int column = 0;

        for (int i = 0; i < nN; i++) {
            //node initialization

            row = i / this.nH;
            column = (i % this.nH);

            nodes[i].x = row * dx;
            nodes[i].y = column * dy;
            nodes[i].T = Data.T0;


            if (row == 0 || column == 0)
            {
                nodes[i].BC = true;
            }
            else if (row == (nB - 1) || column == (nH - 1))
            {
                nodes[i].BC = true;
            }
            else
            {
                nodes[i].BC = false;
            }

        }

        int k = 1;
        for (int i = 0; i < nE; i++)
        {

                // element initialization
                int[] ID =new int[4];

                if (k % nH == 0)
                {
                    //new column
                    k++;
                }

                ID[0] = k;
                ID[1] = (k) + nH;
                ID[2] = (k + nH) + 1;
                ID[3] = (k) + 1;

                k++;

               // System.out.println(ID[0]+" "+ID[1]+" "+ID[2]+" "+ID[3]+" ");

                Node[] nd = new Node[]{nodes[ID[0] - 1], nodes[ID[1] - 1], nodes[ID[2] - 1], nodes[ID[3] - 1]};
                int[] id = new int[] {ID[0], ID[1], ID[2], ID[3]};



                Node[] a = new Node[4];

                for(int l=0;l<4;l++)
                {
                    a[l]=new Node();
                }

    /*
                a[0].x=0;
                a[1].x=0.05;
                a[2].x=0.09;
                a[3].x=0;

                a[0].y=0.;
                a[1].y=0;
                a[2].y=0.04;
                a[3].y=0.02;

                a[0].BC=true;
                a[1].BC=true;
                a[2].BC=true;
                a[3].BC=true;

    */


                elements[i] = new Element(p,nd,id);
       // break;

        }

    }

    void calculateH_global() //agregation
    {
        for (int e = 0; e < nE; e++)
        {
            int[] ID = elements[e].id;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    H_global[ID[i] - 1][ID[j] - 1] += elements[e].H[i][j];
                }
            }
        }
    }

    void calculateHbc_global() //agregation
    {
        for (int e = 0; e < nE; e++)
        {
            int[] ID = elements[e].id;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Hbc_global[ID[i] - 1][ID[j] - 1] += elements[e].Hbc[i][j];
                    H_global[ID[i] - 1][ID[j] - 1] += elements[e].Hbc[i][j];
                }
            }
        }
    }

    void calculateP_global() //agregation
    {
        for (int e = 0; e < nE; e++)
        {
            int[] ID = elements[e].id;

            for (int i = 0; i < 4; i++)
            {
                P_global[ID[i] - 1] += elements[e].P[i];
            }
        }
    }

    void calculateC_global() //agregation
    {
        for (int e = 0; e < nE; e++)
        {
            int[] ID = elements[e].id;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    C_global[ID[i] - 1][ID[j] - 1] += elements[e].C[i][j];
                }
            }
        }
    }

    double[][] calculateC_dt() //[C]/dt
    {
        double[][] C = new double[nN][nN];
        for (int i = 0; i < nN; i++)
        {
            for (int j = 0; j < nN; j++)
            {
                C[i][j] = C_global[i][j] / Data.dt;
            }
        }

        return C;
    }

    double[] calculateP_dt()
    {
        //{P} = {P}+{[C]/dT}*{T0}

        double[] P = new double[nN];
        double[][] C_dt = calculateC_dt();
        for (int i = 0; i < nN; i++)
        {
            P[i] = 0;
            double Cdt = 0;

            for (int j = 0; j < nN; j++)
            {
                Cdt += C_dt[i][j] * nodes[j].T;
            }
            P[i] = P_global[i] + Cdt;
        }


        return P;
    }


    void calculateTemp()
    {
        double[] x = new double[nN];

        //matrix initialization for gauss metod
        double[][] Gauss = new double[nN][nN+1];
        for (int i = 0; i < nN; i++)
        {

        }
       // System.out.println("Iteration "+Data.iteration);

        double[][] C_dt = calculateC_dt();
        double[] P_dt = calculateP_dt();



        for (int i = 0; i < nN; i++)
        {
            for (int j = 0; j < nN; j++)
            {
                Gauss[i][j] = H_global[i][j] + C_dt[i][j];
                //System.out.print(Gauss[i][j]+" ");
            }
            Gauss[i][nN] = P_dt[i];
           // System.out.println(Gauss[i][nN]+" ");
        }

        //calculation gauss elimination
        int l = 1;
        while (l < nN)
        {
            for (int j = 0; j < l; j++)
            {
                if (Gauss[l][j] == 0)
                {
                    //do nothing
                    continue;
                }

                double factor = Gauss[l][j] / Gauss[j][j];
                for (int k = 0; k < (nN + 1); k++)
                {
                    double value_to_substract = 0;

                    if (Gauss[j][k] != 0)
                    {
                        value_to_substract = Gauss[j][k] * factor;
                    }

                    Gauss[l][k] -= value_to_substract;
                }
            }

            l++;
        }



        //calculation of {t1} (claculating variables)
        double[] temp = new double[nN];

        for (int i = nN - 1; i >= 0; i--) {
            temp[i] = 1;

            double value = Gauss[i][nN];
            for (int j = nN - 1; j >= i; j--)
            {
                if (i != j)
                {
                    value -= temp[j] * Gauss[i][j];
                }
                else
                {
                    temp[i] = value / Gauss[i][j];
                }
            }
        }

        //writing results to the {t1}
        for (int i = 0; i < nN; i++) {
            nodes[i].T = temp[i];
        }

        Data.iteration++;
        this.printTemp();

    }


    void printGlobalH()
    {
        System.out.println("Matrix Global H");
        for (int j = 0; j < nN; j++)
        {
            for (int k = 0; k < nN; k++)
            {
                System.out.print(H_global[j][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printGlobalHbc()
    {
        System.out.println("Matrix Global Hbc");
        for (int j = 0; j < nN; j++)
        {
            for (int k = 0; k < nN; k++)
            {
                System.out.print(Hbc_global[j][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printP()
    {
        System.out.println("Matrix Global P");
        for (int j = 0; j < nN; j++)
        {
            System.out.print(P_global[j]+" ");
        }
        System.out.println();
    }

    void printGlobalC()
    {
        System.out.println("Matrix Global C");
        for (int j = 0; j < nN; j++)
        {
            for (int k = 0; k < nN; k++)
            {
                System.out.print(C_global[j][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printElements()
    {
        for(int e = 0; e < nE; e++)
        {
            System.out.println("Elem:"+(e+1)+"\t"+elements[e].id[0]+" "+elements[e].id[1]+" "+elements[e].id[2]+" "+elements[e].id[3]);
        }
        System.out.println();
    }

    void printNodes()
    {
        for (int i = 0; i < nN; i++)
        {
            System.out.println("Nr: "+(i+1)+": "+ nodes[i].x + " " + nodes[i].y + " BC: " + nodes[i].BC);
        }
        System.out.println();
    }

    void printTemp()
    {


        System.out.println("Iteration "+Data.iteration);
        for (int i = 0; i < nN; i++) {
           // System.out.println(nodes[i].T);
        }

        double T_min = Integer.MAX_VALUE;
        double T_max = Integer.MIN_VALUE;


        for (int i = 0; i < nN; i++)
        {
            if (nodes[i].T < T_min)
            {
                T_min = nodes[i].T;
            }
            else
            {
                T_max = nodes[i].T;
            }
        }

        System.out.println("Time: "+Data.iteration*Data.dt+"; T_min:  "+T_min +"; T_max : "+T_max );

    }

    void init()
    {
        for(int i=0;i<this.nN;i++)
        {
            this.nodes[i]=new Node();

        }

        for(int i=0;i<this.nN;i++)
        {
            for(int j=0;j<this.nN;j++)
            {
                this.H_global[i][j]=0;
                this.Hbc_global[i][j]=0;
                this.C_global[i][j]=0;
            }
            this.P_global[i]=0;
        }
    }


}
