import java.util.ArrayList;
import java.util.List;
import static java.lang.Math.*;

public class Wall
{

    PointSchema p;
    double detJ;
    Border border;

    double[][] Hbc = new double[4][4];

    List<double[]> N_wall = new ArrayList<double[]>(); // array of wall shape fucntion
    List<double[]> pc = new ArrayList<double[]>(); // array of wall integration points

    double[] P = new double[4]; //2d array for partial Hbc

    Wall(PointSchema p, Node[] nodes, Border b) //passing 2 element array: start and end point of wall
    {
        this.border = b;
        this.p = p;
        this.Init();


        calculateDetJ(nodes);

        calculateIntegrationPoints();
        calculateN();

        calculateHbc();
        //printHbc();
        calculateP();

    }

    void calculateIntegrationPoints() {
        switch (border)
        {
            case LEFT:
                for (int i = 0; i < p.Psize; i++)
                {
                    pc.add(new double[]{-1,-1 * p.nodes.get(i)}); //set new ksi, eta
                }
                break;
            case DOWN:
                for (int i = 0; i < p.Psize; i++)
                {
                    pc.add(new double[]{p.nodes.get(i),-1}); //set new ksi, eta
                }
                break;
            case RIGHT:
                for (int i = 0; i < p.Psize; i++)
                {
                    pc.add(new double[]{1,-1 * p.nodes.get(i)}); //set new ksi, eta
                }
                break;
            case UP:
                for (int i = 0; i < p.Psize; i++)
                {
                    pc.add(new double[]{p.nodes.get(i),1}); //set new ksi, eta
                }
                break;
            default:
                break;
        }
    }

    void calculateDetJ(Node[] nodes) {
        if (nodes[1].BC == true && nodes[0].BC == true) //jesli istnieje warunek brzegowy na scianie
        {

                double dx=nodes[1].x - nodes[0].x;
                double dy=nodes[1].y - nodes[0].y;
                this.detJ = sqrt(dx*dx + dy*dy)/2;

        }
        else //jeśli niema warunku brzegowego
        {
            this.detJ = 0;
        }
    }

    void calculateN()
    {
        // border initialization with 0
        for (int i = 0; i < p.Psize; i++)
        {
            double[] N_temp= new double[4];
            if (border == Border.LEFT || border == Border.DOWN)
            {
                N_temp[0]= 0.25 * (1 - pc.get(i)[0]) * (1 - pc.get(i)[1]);

            }
            if (border == Border.DOWN || border == Border.RIGHT)
            {
                N_temp[1]= 0.25 * (1 + pc.get(i)[0]) * (1 - pc.get(i)[1]);
            }

            if (border == Border.RIGHT || border == Border.UP)
            {
                N_temp[2]= 0.25 * (1 + pc.get(i)[0]) * (1 + pc.get(i)[1]);
            }
            if (border == Border.UP || border == Border.LEFT)
            {
                N_temp[3]= 0.25 * (1 - pc.get(i)[0]) * (1 + pc.get(i)[1]);
            }
            N_wall.add(N_temp);
        }
    }

    void calculateHbc()
    {
        for (int points = 0; points < p.Psize; points++) //for each point in wall
        {
            double[][] h_temp = new double[4][4];

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    h_temp[i][j] += N_wall.get(points)[i] * N_wall.get(points)[j] * p.weights.get(points) * Data.alpha * detJ;;
                    Hbc[i][j] += h_temp[i][j];
                }
            }
        }
    }


    void calculateP()
    {

        for (int pc = 0; pc < p.Psize; pc++)
        {

            double[] p_temp = new double[]{0,0,0,0};

            for (int i = 0; i < 4; i++)
            {

                p_temp[i] += Data.alpha * p.weights.get(pc) * N_wall.get(pc)[i] * Data.T_ot * detJ;
                P[i] += p_temp[i];
               // System.out.println(P[i]+" "+p_temp[i]+" "+detJ);
            }

        }
        //System.out.println();
    }


    void printHbc()
    {
        System.out.println("Wall" + this.border.toString() );
        System.out.println("Hbc");
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                System.out.print(Hbc[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();

    }


    void Init()
    {
        this.Hbc = new double[][]
                {
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0},
                        {0,0,0,0}
                };
        this.P = new double[] {0,0,0,0};
    }
}

