import java.util.Arrays;
import java.util.List;

import static java.lang.Math.sqrt;

public class P2 extends PointSchema
{

    List<Double> nodesValues = Arrays.asList( -1.0/sqrt(3), 1.0/sqrt(3));
    List<Double> weightValues = Arrays.asList(1.0, 1.0);

    P2()
    {
        this.nodes.addAll(nodesValues);
        this.weights.addAll(weightValues);
        this.Psize = 2;
    }
}