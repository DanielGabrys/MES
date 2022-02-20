import java.util.Arrays;
import java.util.List;

import static java.lang.Math.sqrt;

public class P3 extends PointSchema
{

    List<Double> nodesValues = Arrays.asList(-sqrt(3.0 / 5.0),0.0 , sqrt(3.0 / 5.0));
    List<Double> weightValues = Arrays.asList(5.0/9.0, 8.0/9.0, 5.0/9.0);

    P3()
    {
        this.nodes.addAll(nodesValues);
        this.weights.addAll(weightValues);
        this.Psize = 3;
    }
}