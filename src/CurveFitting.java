import java.io.*;
import java.util.*;

public class CurveFitting {


    static class FastReader {
        BufferedReader br;
        StringTokenizer st;

        public FastReader(String fname) throws FileNotFoundException {
            br = new BufferedReader(new FileReader(fname));
        }

        String next() {
            while (st == null || !st.hasMoreElements()) {
                try {
                    st = new StringTokenizer(br.readLine());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            return st.nextToken();
        }

        int nextInt() {
            return Integer.parseInt(next());
        }

        long nextLong() {
            return Long.parseLong(next());
        }

        double nextDouble() {
            return Double.parseDouble(next());
        }

        String nextLine() throws IOException {
            st = new StringTokenizer(br.readLine());
            String str = "";
            try {
                if (st.hasMoreTokens()) {
                    str = st.nextToken("\n");
                } else {
                    str = br.readLine();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            return str;
        }
    }

    static class Point {
        double xCord;
        double yCord;

        @Override
        public String toString() {
            return "{" +
                    xCord +
                    yCord +
                    '}';
        }

        public Point(double xCord, double yCord) {
            this.xCord = xCord;
            this.yCord = yCord;
        }
    }

    static class Chromosome {
        ArrayList<Double> genes;
        double fitness;

        public Chromosome() {
            this.genes = new ArrayList<>();
            this.fitness = 0.0;
        }

        @Override
        public String toString() {
            return "{" +
                    genes +
                    ",Fitness = " + fitness +
                    '}';
        }
    }

    public static double calcError(Chromosome chromosome, ArrayList<Point> points) {
        double error = 0.0;
        ArrayList<Double> errors = new ArrayList<>();
        for (Point p : points) {
            int power = 0;
            double yCalc = 0.0;
            for (int i = 0; i < chromosome.genes.size(); i++) {
                yCalc += chromosome.genes.get(i) * Math.pow(p.xCord, power);
                power++;
            }
            error = Math.pow((yCalc - p.yCord), 2);
            errors.add(error);
        }
        error = 0.0;
        for (Double e : errors) {
            error += e;
        }
        error *= 1.0 / points.size();
        return error;
    }

    public static int get01Random() {
        return (int) Math.floor(Math.random() * 2);
    }

    public static double getRandomNumber() {
        Random random = new Random();
        return 20.0 * random.nextDouble() - 10.0;
    }

    public static int getRandomNumberSelection(int min, int max) {
        Random random = new Random();
        return random.nextInt(max - min) + min;
    }


    public static double getRandom01Range() {
        Random random = new Random();
        return random.nextDouble();
    }


    static class ChromosomeFitnessComparator implements Comparator<Chromosome> {
        public int compare(Chromosome s1, Chromosome s2) {
            return Double.compare(s1.fitness, s2.fitness);
        }
    }

    public static void calcFitnessSorted(ArrayList<Chromosome> chromosomes, ArrayList<Point> points) {
        for (Chromosome chromosome : chromosomes) {
            chromosome.fitness = calcError(chromosome, points);
        }
        chromosomes.sort(new ChromosomeFitnessComparator());
    }

    public static ArrayList<Chromosome> tournamentSelection(ArrayList<Chromosome> population_chromosomes) {
        ArrayList<Chromosome> tournament = new ArrayList<>();
        for (int i = 0; i < 2; i++) {
            ArrayList<Chromosome> temp = new ArrayList<>();
            int contestants_count = 5;
            while (contestants_count > 0) {
                temp.add(population_chromosomes.get(getRandomNumberSelection(0, population_chromosomes.size())));
                contestants_count--;
            }
            temp.sort(new ChromosomeFitnessComparator());
            tournament.add(temp.get(0));
        }
        return tournament;
    }


    public static void writeToFile(String out){
        try{
            FileWriter writer = new FileWriter("out.txt", true);
            writer.write(out);
            writer.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
    public static void main(String[] args) throws FileNotFoundException {
        File file =new File("out.txt");
        try{
            file.delete();
        }catch (Exception e){
            e.printStackTrace();
        }

        double CROSSOVER_PROBABILITY = 0.6;
        double MUTATION_PROBABILITY = 0.05;
        int MAX_GENERATIONS = 15;
        int T = MAX_GENERATIONS;

        int LB = -10;
        int UB = 10;
        FastReader reader = new FastReader("input.txt");

        int cases = reader.nextInt();

        int index = 0;
        while (cases > 0) {
            int POPULATION_SIZE = 5000;
            ArrayList<Chromosome> chromosomes = new ArrayList<>();
            ArrayList<Chromosome> elite = new ArrayList<>();
            int points = reader.nextInt();
            int degree = reader.nextInt() + 1;
            ArrayList<Point> pointsList = new ArrayList<>();
            while (points > 0) {
                pointsList.add(new Point(reader.nextDouble(), reader.nextDouble()));
                points--;
            }
            //System.out.println(pointsList);

            while (POPULATION_SIZE > 0) {
                int saveDegree = degree;
                Chromosome chromosome = new Chromosome();
                while (saveDegree > 0) {
                    chromosome.genes.add(getRandomNumber());
                    saveDegree--;
                }
                chromosomes.add(chromosome);
                POPULATION_SIZE--;
            }

            //System.out.println(chromosomes);

            int t = 0;
            while (t <= MAX_GENERATIONS) {
                calcFitnessSorted(chromosomes, pointsList);
                int elitePPL = chromosomes.size() / 10;

                while (elitePPL > 0) {
                    elite.add(chromosomes.get(0));
                    chromosomes.remove(chromosomes.get(0));
                    elitePPL--;
                }
                ////////////////Selection/////////////////////

                ArrayList<Chromosome> offsprings = tournamentSelection(chromosomes);
                //  System.out.println(offsprings);


                ///////////////////CrossOver///////////////


                int ITEM_SIZE = chromosomes.get(0).genes.size();
                int cross_over_point1 = getRandomNumberSelection(1, ITEM_SIZE);
                int cross_over_point2 = getRandomNumberSelection(1, ITEM_SIZE);
                // System.out.println(ITEM_SIZE);
                //  System.out.println(cross_over_point1);
                // System.out.println(cross_over_point2);

                while (cross_over_point2 == cross_over_point1) {
                    cross_over_point2 = getRandomNumberSelection(1, ITEM_SIZE);
                }

                int start = Math.min(cross_over_point1, cross_over_point2);
                int end = Math.max(cross_over_point1, cross_over_point2);


                double probability_crossover = getRandom01Range();
                // System.out.println(probability_crossover);

                if (probability_crossover < CROSSOVER_PROBABILITY) {
                    //  System.out.println("Before Cross Over " + offsprings);
                    //   System.out.println("CROSSOVER OCCURRED from point " + start + " To " + end);

                    for (int i = 0; i < ITEM_SIZE; i++) {
                        if (i >= start && i < end) {
                            Double temp = offsprings.get(0).genes.get(i);
                            offsprings.get(0).genes.set(i, offsprings.get(1).genes.get(i));
                            offsprings.get(1).genes.set(i, temp);
                        }
                    }
                    //  System.out.println("After Cross Over " + offsprings);
                } else {
                    //  System.out.println("NO CROSSOVER");
                }

                /////////////Mutation/////////////

                for (int i = 0; i < offsprings.size(); i++) {
                    Chromosome current = offsprings.get(i);
                    // System.out.println("Current Chromosome: " + current);
                    for (int j = 0; j < current.genes.size(); j++) {
                        double Xi = current.genes.get(j);
                        double probability_mutation = getRandom01Range() / 10;
                        if (probability_mutation < MUTATION_PROBABILITY) {
                            double delta_L = Xi - LB;
                            double delta_U = UB - Xi;
                            //   System.out.println("Mutation at gene: " + Xi + " LBi(" + delta_L + ") UBi(" + delta_U + ")");

                            double R1 = getRandom01Range();
                            // System.out.println("R1 Mutation: " + R1);
                            double Y;
                            if (R1 <= 0.5) {
                                Y = delta_L;
                            } else {
                                Y = delta_U;
                            }
                            // System.out.println("Y: " + Y);
                            double r = getRandom01Range();
                            double non_uniformity_deg = 0.5;

                            double p = (1 - ((double) t / T));
                            double p_pow = Math.pow(p, non_uniformity_deg);

                            double delta_t_Y = Y * (1 - Math.pow(r, p_pow));
                            // System.out.println("Delta(t,Y)= " + delta_t_Y);


                            double Xi_new;
                            if (R1 <= 0.5) {
                                Xi_new = Xi - delta_t_Y;
                            } else {
                                Xi_new = Xi + delta_t_Y;
                            }

                            current.genes.set(j, Xi_new);
                            // System.out.println("Xnew= " + Xi_new);
                        }
                    }
                }
                calcFitnessSorted(offsprings, pointsList);
                //System.out.println("Offsprings after Mutation: " + offsprings);

                /////////////Replacement/////////////

                chromosomes.addAll(offsprings);
                chromosomes.addAll(elite);

                chromosomes.sort(new ChromosomeFitnessComparator());

                chromosomes.remove(chromosomes.size() - 1);
                chromosomes.remove(chromosomes.size() - 1);
                //System.out.println("Population: " + chromosomes);
                t++;
            }
            cases--;
            index++;
            System.out.println("\n*********** BEST IS **************** " + chromosomes.get(0));
            StringBuilder out = new StringBuilder();
            out.append("Best for dataset number ").append(index).append("\n");
            out.append("Coefficients:-\n");
            for (int i = 0 ; i < chromosomes.get(0).genes.size();i++){
                out.append(chromosomes.get(0).genes.get(i)).append("\n");
            }
            out.append("Error is ").append(chromosomes.get(0).fitness).append("\n");
            out.append("===============================\n");
            writeToFile(out.toString());

        }
    }
}