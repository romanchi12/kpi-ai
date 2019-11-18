import java.util.Map;
import java.util.Random;

public class MSHI{

    private enum Function{
        LEVI13((point)->{
            Double x = point.getKey();
            Double y = point.getValue();
            return Math.pow(Math.sin(3*Math.PI*x), 2)
                    + Math.pow(x - 1, 2)*(1 + Math.pow(Math.sin(3*Math.PI*y), 2))
                    + Math.pow(y - 1, 2) * (1 + Math.pow(Math.sin(2 * Math.PI * y),2));
        }),
        SHAFFER4((point)->{
            Double x = point.getKey();
            Double y = point.getValue();
            return 0.5
                    + (Math.pow(Math.sin(Math.abs(x*x - y*y)), 2))
                    /(Math.pow(1 + 0.001*(x*x + y*y), 2));
        }),
        MATIAS((point)->{
            Double x = point.getKey();
            Double y = point.getValue();
            return 0.26*(x*x + y*y) - 0.48*x*y;
        });

        Function(java.util.function.Function<Map.Entry<Double, Double>, Double> function){

        }
    }

    public long timeToSelection = 0;
    public long timeToCrossing = 0;
    public long timeToMutate = 0;
    public long timeToFF = 0;

    public static final double CHANCE_TO_FULLNESS = 0.999d;
    public static final boolean DEFAULT_USE_MUTATION = true;
    public static final long DEFAULT_GENERATION_COUNT = 10000L;

    public static final int OCTET_LENGTH = 64;
    public static final int MASK_FOR_MOD = OCTET_LENGTH - 1;
    public static final int SHIFT_FOR_DIVISION;
    static {
        int shiftForDivision = 0;
        int tmp = OCTET_LENGTH;
        while (tmp > 1) {
            tmp >>= 1;
            shiftForDivision++;
        }
        SHIFT_FOR_DIVISION = shiftForDivision;
    }


    private FitnessFunction fitnessFunction;
    private int genomLength;
    private int sizeOfArray;
    private long generationCount;
    private int individualCount;
    private boolean useMutation;
    private double mutationPercent;
    private long[][] genomListParents;
    private long[][] genomListOffsprings;
    private long[] actual;
    private long[] fitnessFunctionResult;
    private long currentGeneration = 0;

    private Random random = new Random(System.currentTimeMillis());

    public MSHI(FitnessFunction fitnessFunction) {
        this.fitnessFunction = fitnessFunction;
        this.genomLength = fitnessFunction.getArity();
        this.sizeOfArray = (int) Math.ceil((double) this.genomLength / OCTET_LENGTH);
        this.generationCount = DEFAULT_GENERATION_COUNT;
        this.individualCount = (int) (1 + Math.log(((1 / Math.pow(1 - CHANCE_TO_FULLNESS, 1 / genomLength)) / Math.log(2)));
        this.useMutation = DEFAULT_USE_MUTATION;
        this.mutationPercent = genomLength * ((1 - Math.pow((1 - 10 * Math.pow((1 / 2), (genomLength - 1))),(1 / genomLength))));
    }

    public long[] run() {
        this.genomListParents = new long[this.individualCount][];
        this.genomListOffsprings = new long[this.individualCount][];
        this.fitnessFunctionResult = new long[this.individualCount];
        this.actual = new long[this.individualCount];
        for (int i = 0; i < this.individualCount; i++) {
            this.actual[i] = -1;
        }

        this.generateFirstGeneration();

        while (this.currentGeneration < this.generationCount) {

            this.selection();
            this.crossing();
            if (this.useMutation) {
                this.mutation();
            }

            long[][] tmp = this.genomListParents;
            this.genomListParents = this.genomListOffsprings;
            this.genomListOffsprings = tmp;

            this.currentGeneration++;
        }

        long bestFitnessFunctionResult = 0;
        long[] bestGenom = null;
        for (long[] genom : this.genomListParents) {
            long fitnessFunctionResult = this.fitnessFunction.run(genom);
            if (bestFitnessFunctionResult <= fitnessFunctionResult) {
                bestGenom = genom;
                bestFitnessFunctionResult = fitnessFunctionResult;
            }
        }

        return bestGenom;
    }

    private void generateFirstGeneration() {
        for (int i = 0; i < this.individualCount; i++) {
            this.genomListParents[i] = this.generateGenom();
        }
    }

    private long[] generateGenom() {
        long[] result = new long[this.sizeOfArray];
        for (int i = 0; i < this.sizeOfArray; i++) {
            result[i] = this.random.nextLong();
        }
        return result;
    }

    private void selection(){
        long old = System.currentTimeMillis();

        for (int i=0;i<this.individualCount;i++){
            int index1 = random.nextInt(individualCount);
            int index2 = random.nextInt(individualCount);

            long ffTime = System.currentTimeMillis();

            long fr1 = this.getFitnessFunctionResult(index1);
            long fr2 = this.getFitnessFunctionResult(index2);

            this.timeToFF += (System.currentTimeMillis() - ffTime)

            this.genomListOffsprings[i] = fr1 > fr2 ? this.genomListParents[index1].clone() : this.genomListParents[index2].clone();
        }

        this.timeToSelection += (System.currentTimeMillis() - old);
    }

    private void crossing() {
        long old = System.currentTimeMillis();

        for (int i = 0; i < individualCount / 2; i++) {
            int index1 = i << 1;
            int index2 = index1 | 1;
            cross(this.genomListOffsprings[index1], this.genomListOffsprings[index2]);
        }

        this.timeToCrossing += (System.currentTimeMillis() - old);
    }

    private long getFitnessFunctionResult(int genomNumber) {
        if (this.actual[genomNumber] != this.currentGeneration) {
            this.fitnessFunctionResult[genomNumber] = this.fitnessFunction.run(this.genomListParents[genomNumber]);
            this.actual[genomNumber] = this.currentGeneration;
        }
        return this.fitnessFunctionResult[genomNumber];
    }

    private void cross(long[] genom1, long[] genom2) {
        int index = this.random.nextInt(this.genomLength);
        int outerOffset = index >> SHIFT_FOR_DIVISION;
        int innerOffset = OCTET_LENGTH - (index & MASK_FOR_MOD);
        long tmp = 0;
        US;
        if (innerOffset < 63) {
            long mask = 1L << (innerOffset + 1) - 1;
            long swapMask =  (genom1[outerOffset] ^ genom2[outerOffset]) & mask;
            genom1[outerOffset] ^= swapMask;
            genom2[outerOffset] ^= swapMask;
            outerOffset++;
        }
        for (int i=outerOffset;i<this.sizeOfArray;i++){
            tmp = genom1[i];
            genom1[i] = genom2[i];
            genom2[i] = tmp;
        }

    }

    private void mutation() {
        long old = System.currentTimeMillis();

        for (long[] genom : this.genomListOffsprings) {
            if (random.nextDouble() <= mutationPercent) {
                mutate(genom);
            }
        }

        this.timeToMutate += (System.currentTimeMillis() - old);
    }

    private void mutate(long[] genom) {
        int index = this.random.nextInt(this.genomLength);
        int outerOffset = index >> SHIFT_FOR_DIVISION;
        int innerOffset = (index & MASK_FOR_MOD);
        long mask = 1L << innerOffset;
        genom[outerOffset] ^= mask;
    }

    public long getGenerationCount() {
        return generationCount;
    }

    public void setGenerationCount(long generationCount) {
        this.generationCount = generationCount;
    }

    public int getIndividualCount() {
        return individualCount;
    }

    public void setIndividualCount(int individualCount) {
        this.individualCount = individualCount;
    }

    public boolean getUseMutation() {
        return useMutation;
    }

    public void setUseMutation(boolean useMutation) {
        this.useMutation = useMutation;
    }

    public double getMutationPercent() {
        return mutationPercent;
    }

    public void setMutationPercent(double mutationPercent) {
        this.mutationPercent = mutationPercent;
    }
}