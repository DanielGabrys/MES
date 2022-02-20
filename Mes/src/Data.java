import java.io.*;
import java.util.Scanner;

public class Data
{
    static double T0; //initial temperature
    static  double dt; //simulation step time
    static  double t; //simulation time
    static double T_ot; //ambient temperature
    static  double alpha; //alfa


    static   double H; //height
    static   double B; //width
    static  int nH; //nodes H
    static  int nB; //nodes B

    static double c; //specific heat
    static double k; //conductivity
    static double ro; //density

    static int iteration=0; //current iteration


    public void dataRead(String filename)
    {

        int variables =12;
        int counter=0;
        double[] data= new double[variables];
        try
        {
            Scanner br = new Scanner(new FileReader("src/"+filename));
            while (counter<variables)
            {
                //System.out.println(br.nextLine()+" "+counter);
                data[counter]=Double.parseDouble(br.nextLine());

                counter++;
            }
            br.close();
        }
        catch (FileNotFoundException e)
        {
            System.err.println("File not found");
        }

        for(int i=0;i<variables;i++)
        {
            switch(i)
            {
                case 0:
                {
                    this.T0=data[0];
                }
                break;
                case 1:
                {

                    this.dt=data[1];
                }
                break;
                case 2:
                {

                    this.t=data[2];
                }
                break;
                case 3:
                {

                    this.T_ot=data[3];
                }
                break;
                case 4:
                {

                    this.alpha=data[4];
                }
                break;
                case 5:
                {

                    this.H=data[5];
                }
                break;
                case 6:
                {

                    this.B=data[6];
                }
                break;
                case 7:
                {

                    this.nH=(int) data[7];
                }
                break;
                case 8:
                {
                    this.nB=(int) data[7];
                }
                break;
                case 9:
                {
                    this.c=data[9];
                }
                break;
                case 10:
                {
                    this.k=data[10];
                }
                break;
                case 11:
                {
                    this.ro=data[11];
                }
                break;
            }
        }
}

}
