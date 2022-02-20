public class Main
{
    public static void main(String[] args)
    {
       Data load_file = new Data();
       //load_file.dataRead("data.txt");
       load_file.dataRead("data2.txt");

        /*
        System.out.println(Data.T0);
        System.out.println(Data.dt);
        System.out.println(Data.t);
        System.out.println(Data.T_ot);
        System.out.println(Data.alpha);
        System.out.println(Data.H);
        System.out.println(Data.B);
        System.out.println(Data.nH);
        System.out.println(Data.nB);
        System.out.println(Data.c);
        System.out.println(Data.k);
        System.out.println(Data.ro);
        */



       PointSchema P = new PointSchema();
       P2 p2 = new P2();
       P3 p3 = new P3();
       Grid g = new Grid(p3);



       g.printNodes();
       g.printElements();

       g.calculateH_global();
       g.calculateHbc_global();
       g.calculateP_global();
       g.calculateC_global();

       for(int i=0;i<Data.t/Data.dt;i++)
       {
           g.calculateTemp();
       }





    }


}
