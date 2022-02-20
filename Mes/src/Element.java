import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Element
{
    int size;
    int [] id;
    Node [] nodes;
    Jakobian[] J;
    PointSchema Pt;

    List <double[]> ShapeValues= new ArrayList<double[]>(); //2d array for [integration_poits][N]
    List <double[]> ShapeValuesDKsi= new ArrayList<double[]>(); //2d array for [integration_poits][Ndksi]
    List <double[]> ShapeValuesDEta= new ArrayList<double[]>(); //2d array for [integration_poits][NDeta]

    List <double[]> dNx= new ArrayList<double[]>(); //2d array list for dNX
    List <double[]> dNy= new ArrayList<double[]>(); //2d array for dNY

    double[][] H= new double[4][4]; // array 2D H local
    double[][] Hbc= new double[4][4]; // array 2D Hbc local
    Wall[] walls;

    double[] P= new double[4]; //2d array for partial Hbc
    double[][] C= new double[4][4]; // array 2D C local


    Element(PointSchema p,Node[] nodes,int[] id )
    {
        this.nodes = nodes;
        this.size = p.Psize;
        this.id = id;
        J =new Jakobian[size*size];
        this.Pt=p;

        this.Init();
        
        //temp arrays
        double[] tem_N_arr = new double[4]; //temporary array for N shape function values
        double[] tem_Nde_arr = new double[4]; //temporary array for pochodna Eta shape function values
        double[] tem_Ndk_arr = new double[4]; //temporary array for N shape Ksi function values

        int c=0;

        // calculate N,Ndksi, Ndeta
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {

                double ksi = p.nodes.get(j);
                double eta = p.nodes.get(i);
                // ksi = -0.23;
                // eta = 0.55;

                // wartości funkcji kształtu dla c-tego punktu całkowania
                tem_N_arr[0] = 0.25 * (1 - ksi) * (1 - eta);
                tem_N_arr[1] = 0.25 * (1 + ksi) * (1 - eta);
                tem_N_arr[2] = 0.25 * (1 + ksi) * (1 + eta);
                tem_N_arr[3] = 0.25 * (1 - ksi) * (1 + eta);


                ShapeValues.set(c,new double[] {tem_N_arr[0],tem_N_arr[1],tem_N_arr[02],tem_N_arr[3], });

                /*
                for(int k=0;k<4;k++)
                     System.out.print(ShapeValues.get(c)[k]+" ");
                System.out.println();
                */

                // wartości pochodnych funkcji kształtu po ksi dla c-tego punktu całkowania
                tem_Ndk_arr[0] = -0.25 * (1 - eta);
                tem_Ndk_arr[1] = 0.25 * (1 - eta);
                tem_Ndk_arr[2] = 0.25 * (1 + eta);
                tem_Ndk_arr[3] = -0.25 * (1 + eta);


                ShapeValuesDKsi.set(c,new double[]{tem_Ndk_arr[0],tem_Ndk_arr[1],tem_Ndk_arr[2],tem_Ndk_arr[3]});

                // wartości pochodnych funkcji kształtu po eta dla c-tego punktu całkowania
                tem_Nde_arr[0] = -0.25 * (1 - ksi);
                tem_Nde_arr[1] = -0.25 * (1 + ksi);
                tem_Nde_arr[2] = 0.25 * (1 + ksi);
                tem_Nde_arr[3] = 0.25 * (1 - ksi);

                ShapeValuesDEta.set(c,new double[]{tem_Nde_arr[0],tem_Nde_arr[1],tem_Nde_arr[2],tem_Nde_arr[3]});

                c++;
            }
        }


        // calculate Jakobian for each integration point
        for(int i=0;i<size*size;i++)
        {
            calculateJ(i);
            calculateDxDy();
        }

        wallInit();
        calculateH(p);
        calculateHbc();
        calculateP();
        calculateC(p);


        //printDNksi();
        //printDNeta();
        //printInvJacobian(0);
        //printDNx();
        //printDNy();
        //printH();
        //printC();




    }

    void calculateJ(int i)
    {

        Jakobian jak = new Jakobian();

        //1 row of J matrix
        double dXdKsi = 0.0;
        double dYdKsi = 0.0;

        //2 row of  J matrix
        double dXdEta = 0.0;
        double dYdEta = 0.0;


        double[] X = {nodes[0].x, nodes[1].x, nodes[2].x, nodes[3].x}; // x1, x2 ,x3,x4
        double[] Y= {nodes[0].y, nodes[1].y, nodes[2].y, nodes[3].y};

        double x=0;
        //interpolation
        for (int j = 0; j < 4; j++)
        {
            dXdKsi += this.ShapeValuesDKsi.get(i)[j] * X[j];
            dXdEta += this.ShapeValuesDEta.get(i)[j] * X[j];
            dYdKsi += this.ShapeValuesDKsi.get(i)[j] * Y[j];
            dYdEta += this.ShapeValuesDEta.get(i)[j] * Y[j];

            x += this.ShapeValues.get(i)[j] * X[j];
            //System.out.print(this.ShapeValues.get(i)[j]+"*"+X[j]+"+");
        }
        //System.out.println();
        //System.out.println(x);


        jak.matrixJ[0][0] = dXdKsi;
        jak.matrixJ[0][1] = dYdKsi;
        jak.matrixJ[1][0] = dXdEta;
        jak.matrixJ[1][1] = dYdEta;

        //detJ
        jak.detJ = jak.matrixJ[0][0] * jak.matrixJ[1][1] - jak.matrixJ[0][1] * jak.matrixJ[1][0];


        jak.inv_m[0][0] = 1.0 / jak.detJ * jak.matrixJ[1][1];
        jak.inv_m[0][1] = -1.0 / jak.detJ * jak.matrixJ[0][1];
        jak.inv_m[1][0] = -1.0 / jak.detJ * jak.matrixJ[1][0];
        jak.inv_m[1][1] = 1.0 / jak.detJ * jak.matrixJ[0][0];

        J[i] = jak;
        //System.out.println(J[i].detJ);

    };

    void calculateDxDy()
    {
        //for each integration point

        for (int i = 0; i < size * size; i++)
        {
            double[] dNdX=new double[4];
            double[] dNdY=new double[4];

            //for every shape function
            for (int j = 0; j < 4; j++)
            {
                dNdX[j] = J[i].inv_m[0][0] * ShapeValuesDKsi.get(i)[j] + J[i].inv_m[0][1] * ShapeValuesDEta.get(i)[j];
                dNdY[j] = J[i].inv_m[1][0] * ShapeValuesDKsi.get(i)[j] + J[i].inv_m[1][1] * ShapeValuesDEta.get(i)[j];

            }
            dNx.set(i,dNdX);
            dNy.set(i,dNdY);

        }
    }

    void calculateH(PointSchema n)
    {
    //for each integration point
    for(int i = 0; i < size *  size; i++)
    {
        int row = i / size;
        int column = i % size;

        double[][] h_temp = new double[4][4];

        //calculate Hp for point
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                //dx
                h_temp[j][k] = dNx.get(i)[j] * dNx.get(i)[k];

                //dy
                h_temp[j][k] += dNy.get(i)[j] * dNy.get(i)[k];

                h_temp[j][k] = h_temp[j][k] * Data.k * J[i].detJ * n.weights.get(row) * n.weights.get(column); //weights

                // sum right Hpc matric fields to local matrix H
                H[j][k] += h_temp[j][k];
            }
        }
        //printHBPC(h_temp);
    }

}

    void calculateHbc() //sum wall hbc matrixes
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++) {
                    Hbc[j][i] += walls[k].Hbc[j][i];
                }
            }
    }

    }

    void calculateP()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                P[i] += walls[j].P[i];

            }
        }
    }


    void calculateC(PointSchema n)
    {
        for (int i = 0; i < size * size; i++)
        {
            int row = i / size; //which row
            int column = i % size; //which column

            double[][] C_temp = new double[4][4];

            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    C_temp[j][k] += ShapeValues.get(i)[j] * ShapeValues.get(i)[k];
                    C_temp[j][k] = C_temp[j][k] * Data.c * Data.ro * J[i].detJ * n.weights.get(row) * n.weights.get(column);

                    C[j][k] += C_temp[j][k];
                }
            }
        }
    }


    void printH()
    {

        System.out.println("H matrix: ");

        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                System.out.print(H[j][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printDNx()
    {

        System.out.println("dNdX");
        for (int i = 0; i < size * size; i++)
        {
            for (int j = 0; j < 4; j++) {

                System.out.print(dNx.get(i)[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
}

    void printDNy()
    {
        System.out.println("dNdY");
        for (int i = 0; i < size * size; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(dNy.get(i)[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printDNksi()
    {
        System.out.println("dNdksi");
        for (int i = 0; i < size * size; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                System.out.print(ShapeValuesDKsi.get(i)[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printDNeta()
    {
        System.out.println("dNdEta");
        for (int i = 0; i < size * size; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(ShapeValuesDEta.get(i)[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printHBPC(double [][] hp)
    {
        System.out.println("HBPC");
        for (int i = 0; i < size * size; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                System.out.print(hp[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void printInvJacobian(int point)
    {
    System.out.println("INv Jakobian");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            System.out.print(J[point].inv_m[i][j] + " ");
        }
        System.out.println();
    }
    System.out.println();
}

    void printC()
    {
    System.out.println("Matrix C");
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                System.out.print(C[j][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    void Init()
    {
        double[] temp = new double[] {0,0,0,0};
        double[][] temp2 = new double[4][4];
        for(int i=0;i<size*size;i++)
        {
            if(i<4)
            {
                temp2[i]=temp;
            }
                this.ShapeValues.add(temp);
                this.ShapeValuesDEta.add(temp);
                this.ShapeValuesDKsi.add(temp);
                this.dNx.add(temp);
                this.dNy.add(temp);
                this.J[i] = new Jakobian();

        }
        this.H = new double[][]
                {
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0}
                };
        this.C= new double[][]
                {
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0}
                };

    }

    void wallInit()
    {

        walls = new Wall[4];
        walls[0] = new Wall(this.Pt, new Node[]{ nodes[0], nodes[1] }, Border.DOWN);
        walls[1] = new Wall(this.Pt, new Node[]{ nodes[1], nodes[2] }, Border.RIGHT);
        walls[2] = new Wall(this.Pt, new Node[]{ nodes[3], nodes[2] }, Border.UP);
        walls[3] = new Wall(this.Pt, new Node[]{ nodes[0], nodes[3] }, Border.LEFT);
    }
}
