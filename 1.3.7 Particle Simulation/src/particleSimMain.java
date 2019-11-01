import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;

// figuring out if an enzyme will drift to high or low substrate concentration

/**
 * <body1>This is a particle simulation program that attempts to simulate the
 * movements of 2 particles connected by a spring in 3D space.</Body1><br>
 * 
 *
 * 
 * @author Nathan Pollock
 * @version V1.3.7
 *          <h1>Patch notes since V1.3.1</h1>- Fixed basic physics
 *          <h2>Planned Improvements</h2>- test random forces<br>
 *          - clean and comment code
 * @since 6-26-2019
 */

//energy remains constant/ low oscillation
//added friction
//added random force

//added 
public class particleSimMain {
	// Variables:
	// public static BufferStrategy bs;
	// public static Graphics g;
	// doubles:
	public static double gausMult, rms, startTime, currLegnth, p1AvgVel, p2AvgVel;
	public static double KE, PE;
	// ArrayLists:
	public static ArrayList<Particle> particles;
	public static ArrayList<Position> midpoints;
	public static ArrayList<Double> energies;
	public static ArrayList<Double> PEenergy;
	public static ArrayList<Double> KEenergy;
	public static ArrayList<Double> bondLegnths;
	public static ArrayList<Double> midpointXs;
	// Particles:
	public static Particle particle1, particle2;
	// booleans:
	public static boolean firstRun = true;
	// ints:
	public static int currStep = 0;

	// ======================= Main Settings ========================
	public static double bondLegnth = 1.54 * Math.pow(10, -10); // in m
	public static double initialBondLegnth = 1.54 * Math.pow(10, -10);
	public static double springConstant = 2.92; // kg/sec^2
	public static double frictionCoef = 6 * Math.PI * 8.9 * Math.pow(10, -4) * 1.5 * Math.pow(10, -10); // kg/s
	public static double particleMass = (12.01 / (6.022 * Math.pow(10, 23)) / 1000); // in kg
	public static double substrateGradientCoeff = Math.pow(10, -11); // particles/cm^2
	public static boolean constant = true;
	public static double xOffset = 0;

	// box boundaries (+ or -)
	public static double xBound = Math.pow(10, -9);
	public static double yBound = Math.pow(10, -9);
	public static double zBound = Math.pow(10, -9);

	//public static double targetTemp = 300;// k//Math.pow(10, 52); // kelvin
	public static double TMPgauMult = Math.pow(10, -8);
	public static double TMPMax = 400;
	public static double TMPLo = 100;

	public static boolean showAvgXYZ = false;
	public static boolean showMidpointXAllConsole = false;
	public static boolean showDistanceData = false;
	public static boolean showBondLegnthData = false;
	public static boolean printVelocities = false;
	public static boolean showEnergyDataAllConsole = false;
	public static boolean showPEDataAllConsole = false;
	public static boolean showKEDataAllConsole = false;
	public static boolean showKEDataAverages = false;
	public static boolean saveEnergyDataAll = false;
	public static boolean showEnergyDataSimple = false;
	public static boolean printForces = false;
	public static boolean showEnergyDataGraph = false;
	public static boolean showEnergyAveragesGraph = false;
	public static boolean showMidpointXGraph = false;
	public static boolean showP1XGraph = false;
	public static boolean showAvgPositionX = true;

	public static boolean randomForce = true;
	public static boolean springForce = true;
	public static boolean friction = true;
	public static boolean substrateGradient = true;

	public static boolean setInitialVelocity = false;
	public static boolean setInitialForce = false;
	public static boolean boundaries = true;

	public static boolean calculateGausNew = false;

	public static boolean saveData = false;
	public static boolean openVMD = false;
	public static boolean copyToClipboard = false;

	public static int timeSteps = 1000000;
	public static double scale = Math.pow(10, -15); // 1 calculation per femtosecond is standard

	public static double desentPerc = .01;

	// sets an initial force on particle1 for the first time step
	public static Position initialVelocity = new Position(100, 0, 0);

	public static Direction initialForceDir = new Direction(1, 0, 0);
	public static Double initialForceMag = (double) 1; // newton

	public static String path = "C:\\Users\\npoll\\AppData\\Roaming\\Microsoft\\Windows\\Start Menu\\Programs\\University of Illinois\\VMD\\particleSimulationOutput.xyz";
	public static String dataOutPath = "C:\\Users\\npoll\\Documents\\particleData.txt";
	public static String dataOutPath2 = "C:\\Users\\npoll\\Documents\\particleData2.txt";
	public static String dataOutPath3 = "C:\\Users\\npoll\\Documents\\particleData3.txt";
	public static String energyDataOutPath = "C:\\Users\\npoll\\Documents\\energyDataOutput.txt";

	// ==================== End of Main Settings =====================
	public static void main(String[] args) {
		init();
		mainPhase();
		printOutsPhase();
		if (saveData) {
			save();
		}
		if (copyToClipboard) {
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(path), null);
		}
		if (openVMD) {
			try {
				Runtime.getRuntime().exec("cmd /c start \"\" C:\\Users\\npoll\\Desktop\\vmd.bat");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		graphingPhase();
	}

	public static void graphingPhase() {
		if (showMidpointXGraph) {

			for (int g = 0; g < midpoints.size(); g++) {
				midpointXs.add(midpoints.get(g).x);
			}
			new Graph("Midpoint X vs Time", 800, 400, midpointXs);
		}
		if (showEnergyDataGraph) {
			new Graph("Energy Vs Time", 800, 400, energies);
		}
		if (showP1XGraph) {
			ArrayList<Double> xp1 = new ArrayList<Double>();
			for (int i = 0; i < particle1.positions.size(); i++) {
				xp1.add(particle1.positions.get(i).x);
			}
			new Graph("Particle 1 X Vs Time", 800, 400, xp1);
		}

		if (showEnergyAveragesGraph) {
			int binSize = 100;
			double binSum = 0;
			ArrayList<Double> averages = new ArrayList<Double>();
			for (int iter = 0; iter < timeSteps / binSize; iter++) {
				for (int y = 0; y < binSize * iter; y++) {
					binSum += KEenergy.get(y);

				}
				averages.add(binSum / (binSize * iter));
				binSum = 0;
			}
			new Graph("Average Energy Vs Time", 800, 400, averages);
		}
		if(showAvgPositionX) {
			double sum = 0;
			for(int i = 0; i<midpoints.size(); i++) {
				sum+=midpoints.get(i).x;
			}
			System.out.println(sum/midpoints.size());
		}
	}

	public static void printOutsPhase() {

		if (showAvgXYZ) {
			System.out.println("============= XYZ Averages ============");
			double xSum = 0;
			double ySum = 0;
			double zSum = 0;
			String content1 = "";
			for (int i = 0; i < timeSteps - 1; i++) {
				content1 += (midpoints.get(i + 1).x - midpoints.get(i).x) + "\n";
				xSum += (midpoints.get(i + 1).x - midpoints.get(i).x);
				ySum += (midpoints.get(i + 1).y - midpoints.get(i).y);
				zSum += (midpoints.get(i + 1).z - midpoints.get(i).z);

			}

			System.out.println("X Average: " + xSum / (timeSteps - 1));
			System.out.println("Y Average: " + ySum / (timeSteps - 1));
			System.out.println("Z Average: " + zSum / (timeSteps - 1));
			System.out.println("Saving X data...");
			writeToFile(dataOutPath2, content1);
			System.out.println("Saved");
		}
		if (showMidpointXAllConsole) {
			for (int g = 0; g < midpoints.size(); g++) {
				System.out.println(midpoints.get(g).x);
				midpointXs.add(midpoints.get(g).x);
			}
		}
		if (showDistanceData) {
			System.out.println("=========== Distance Data ===========");
			String content = "";
			for (int k = 0; k < timeSteps - 1; k++) {
				content += distance(midpoints.get(k), midpoints.get(k + 1)) + "\n";

				System.out.println(distance(midpoints.get(k), midpoints.get(k + 1)));
			}
			System.out.println("Saving distance data...");
			writeToFile(dataOutPath, content);
			System.out.println("Saved");
		}
		if (showEnergyDataAllConsole) {
			for (int p = 0; p < energies.size(); p++) {
				System.out.println(energies.get(p));
			}
		}
		if (showPEDataAllConsole) {
			for (int q = 0; q < PEenergy.size(); q++) {
				System.out.println(PEenergy.get(q));
			}
		}
		if (showKEDataAllConsole) {
			for (int w = 0; w < energies.size(); w++) {
				System.out.println(energies.get(w));
			}
		}
		if (saveEnergyDataAll) {
			for (int o = 0; o < energies.size(); o++) {
				energies.get(o);
			}

		}
		if (showEnergyDataSimple) {
			System.out.println("============== Simple Energy Data ===============");
			double energyAvg = alistAvg(energies);
			extrema min = min(energies);
			extrema max = max(energies);
			double range = max.magnitude - min.magnitude;
			double percent = 100 * (range / energyAvg);
			System.out.println("Average Energy: " + energyAvg + "\nMax Energy: " + max.magnitude + " on frame "
					+ max.timeStep + "\nMin Energy: " + min.magnitude + " on frame " + min.timeStep);
			System.out.println("Range: " + range + "\nPercent Change: " + percent + "%");
		}
		if (showBondLegnthData) {
			for (int d = 0; d < bondLegnths.size(); d++) {
				System.out.println(bondLegnths.get(d));
			}
		}
		if (showKEDataAverages) {
			int binSize = 100;
			double binSum = 0;
			ArrayList<Double> averages = new ArrayList<Double>();
			for (int iter = 0; iter < timeSteps / binSize; iter++) {
				for (int y = 0; y < binSize * iter; y++) {
					binSum += KEenergy.get(y);

				}
				averages.add(binSum / (binSize * iter));
				System.out.println(averages.get(averages.size() - 1));
				binSum = 0;
			}

		}

	}

	private static void init() {
		startTime = System.currentTimeMillis();
		particles = new ArrayList<Particle>();
		midpoints = new ArrayList<Position>();
		energies = new ArrayList<Double>();
		PEenergy = new ArrayList<Double>();
		KEenergy = new ArrayList<Double>();
		bondLegnths = new ArrayList<Double>();
		midpointXs = new ArrayList<Double>();
		// instantiate particles
		particle1 = new Particle(new Position((-initialBondLegnth / 2) + xOffset, 0, 0), "C");
		particle2 = new Particle(new Position((initialBondLegnth / 2) + xOffset, 0, 0), "C");
		if (setInitialVelocity) {
			particle1.velocity = new Position(initialVelocity.x, initialVelocity.y, initialVelocity.z);
		} else {
			particle1.velocity = new Position(0, 0, 0);
		}
		particle2.velocity = new Position(0, 0, 0);
		// add particles to particles array list
		particles.add(particle1);
		particles.add(particle2);
		// get Target 1-Dimensional velocity - unused
		// rms = (double) Math.sqrt((double) (3 * 8.31 * targetTemp / particleMass));
		gausMult = (getGausMult(rms, calculateGausNew));
	}

	private static void mainPhase() {
		// here we generate the x,y,z positions of each particle during each time stamp
		// to be exported to VMD in an .xyz file format
		while (currStep < timeSteps) {
			// main loop start
			updateParticleForces();
			updateParticleVelocities();
			updateParticlePositions();
			recordPositions();
			recordEnergies();
			currStep++;
			firstRun = false;
		} // end of main loop

	}

	private static void recordEnergies() {

		PE = .5 * springConstant * Math.pow(Math.abs(currLegnth - bondLegnth), 2);
		KE = (.5 * particleMass * Math.pow(p1AvgVel, 2)) + (.5 * particleMass * Math.pow(p2AvgVel, 2));

		PEenergy.add(PE);
		KEenergy.add((KE));

		energies.add((KEenergy.get(KEenergy.size() - 1) + PEenergy.get(PEenergy.size() - 1)));
	}

	public static void updateParticlePositions() {

		particle1.moveVec(mult(particle1.velocity, scale));
		particle2.moveVec(mult(particle2.velocity, scale));

		if (boundaries) {

			// above bound p1
			if (particle1.pos.x > xBound) {
				particle1.setPos(new Position(xBound, particle1.pos.y, particle1.pos.z));
			}
			if (particle1.pos.y > yBound) {
				particle1.setPos(new Position(particle1.pos.x, yBound, particle1.pos.z));
			}
			if (particle1.pos.z > zBound) {
				particle1.setPos(new Position(particle1.pos.x, particle1.pos.y, zBound));
			}
			// below bound p1
			if (particle1.pos.x < -xBound) {
				particle1.setPos(new Position(-xBound, particle1.pos.y, particle1.pos.z));
			}
			if (particle1.pos.y < -yBound) {
				particle1.setPos(new Position(particle1.pos.x, -yBound, particle1.pos.z));
			}
			if (particle1.pos.z < -zBound) {
				particle1.setPos(new Position(particle1.pos.x, particle1.pos.y, -zBound));
			}

			// above bound p2
			if (particle2.pos.x > xBound) {
				particle2.setPos(new Position(xBound, particle2.pos.y, particle2.pos.z));
			}
			if (particle2.pos.y > yBound) {
				particle2.setPos(new Position(particle2.pos.x, yBound, particle2.pos.z));
			}
			if (particle2.pos.z > zBound) {
				particle2.setPos(new Position(particle2.pos.x, particle2.pos.y, zBound));
			}
			// below bound p2
			if (particle2.pos.x < -xBound) {
				particle2.setPos(new Position(-xBound, particle2.pos.y, particle2.pos.z));
			}
			if (particle2.pos.y < -yBound) {
				particle2.setPos(new Position(particle2.pos.x, -yBound, particle2.pos.z));
			}
			if (particle2.pos.z < -zBound) {
				particle2.setPos(new Position(particle2.pos.x, particle2.pos.y, -zBound));
			}

		}

		if (printForces) {
			System.out.println("P1 Pos: (" + (float) particle1.pos.x + "," + (float) particle1.pos.y + ","
					+ (float) particle1.pos.z + ")");
			System.out.println("P2 Pos: (" + (float) particle2.pos.x + "," + (float) particle2.pos.y + ","
					+ (float) particle2.pos.z + ")");
		}

	}

	public static void updateParticleVelocities() {
		// vec sum takes care of scaling
		double p1PrevVel = vecToMag(particle1.velocity);
		double p2PrevVel = vecToMag(particle2.velocity);

		particle1.setVelocity(add(particle1.velocity, mult(divide(vecSum(particle1.forces), particleMass), scale)));
		particle2.setVelocity(add(particle2.velocity, mult(divide(vecSum(particle2.forces), particleMass), scale)));

		p1AvgVel = (p1PrevVel + vecToMag(particle1.velocity)) / 2;
		p2AvgVel = (p2PrevVel + vecToMag(particle2.velocity)) / 2;

		if (printVelocities) {
			System.out.println("Frame: " + currStep + "\nParticle1 Velocity (m/sec): (" + particle1.velocity.x + ","
					+ particle1.velocity.y + "," + particle1.velocity.z + "\nParticle2 Velocity (m/sec): ("
					+ particle2.velocity.x + "," + particle2.velocity.y + "," + particle2.velocity.z + ")\n");
		}

	}

	public static void updateParticleForces() {

		for (int i = 0; i < particles.size(); i++) {
			// clear all forces
			particles.get(i).forces.removeAll(particles.get(i).forces);
		}
		// current bond length of atoms

		if (firstRun && setInitialForce) {
			particle1.forces.add(new Vector(initialForceDir, initialForceMag, "Initial Force"));
		}

		// spring force calculation
		currLegnth = distance(particle1.pos, particle2.pos);
		bondLegnths.add(currLegnth);
		Position midpoint = midpoint(particle1.pos, particle2.pos);

		Direction p1Direction = sub(particle2.pos, particle1.pos).toDir();
		double p1Mag = (double) ((springConstant * (currLegnth - bondLegnth)));

		Direction p2Direction = sub(particle1.pos, particle2.pos).toDir();
		double p2Mag = (double) ((springConstant * (currLegnth - bondLegnth)));

		Vector p1Force = (new Vector(p1Direction, p1Mag, "Spring"));
		Vector p2Force = new Vector(p2Direction, p2Mag, "Spring");

		if (springForce) {
			particle1.forces.add(p1Force);
			particle2.forces.add(p2Force);
		}
		// end of spring force calculation

		// friction force calculation

		Vector p1Drag = new Vector(mult(particle1.velocity, -1).toDir(), vecToMag(particle1.velocity) * frictionCoef,
				"Friction");
		Vector p2Drag = new Vector(mult(particle2.velocity, -1).toDir(), vecToMag(particle2.velocity) * frictionCoef,
				"Friction");

		if (friction) {
			particle1.forces.add(p1Drag);
			particle2.forces.add(p2Drag);
		}
		// end of friction force calculation

		Vector p1Random, p2Random;
		if (!substrateGradient) {
			p1Random = new Vector(new Direction(getRand(), getRand(), getRand()), gau(gausMult), "Random Force");
			p2Random = new Vector(new Direction(getRand(), getRand(), getRand()), gau(gausMult), "Random Force");
		} else {
			p1Random = new Vector(new Direction(getRand(), getRand(), getRand()),
					getMult(midpoint.x), "Random Force");
			p2Random = new Vector(new Direction(getRand(), getRand(), getRand()),
					getMult(midpoint.x), "Random Force");
		}
		if (randomForce) {
			particle1.forces.add(p1Random);
			particle2.forces.add(p2Random);
		}
		// end of random force calculation

		midpoint = midpoint(particle1.pos, particle2.pos);
		midpoints.add(midpoint);

		if (printForces) {
			System.out.println("Frame: " + currStep);
			for (int i = 0; i < particle1.forces.size(); i++) {
				System.out.println(
						"P1 Force " + i + " " + particle1.forces.get(i).tag + ": " + particle1.forces.get(i).magnitude);
			}
			for (int k = 0; k < particle2.forces.size(); k++) {
				System.out.println(
						"P2 Force " + k + " " + particle2.forces.get(k).tag + ": " + particle2.forces.get(k).magnitude);
			}
		}

	}

	public static void save() {
		System.out.println();
		System.out.println("Saving main data...");
		writeToFile(path, getText());
		System.out.println("Saved to " + path);
		System.out.println("Completed " + timeSteps + " timeSteps of " + particles.size() + " particles in "
				+ (System.currentTimeMillis() - startTime) + "ms");
	}

	public static String getText() {
		int inc = 1;
		int n = 1000;
		String content = "";
		if (timeSteps > n && timeSteps % n == 0) {
			inc = timeSteps / n;
		}
		System.out.println("Saving every " + inc + " timestep");
		for (int i = 0; i < timeSteps; i += inc) {
			content += particles.size() + "\n" + i;
			for (int j = 0; j < particles.size(); j++) {
				content += "\n" + particles.get(j).type + " " + particles.get(j).positions.get(i).x * Math.pow(10, 10)
						+ " " + particles.get(j).positions.get(i).y * Math.pow(10, 10) + " "
						+ particles.get(j).positions.get(i).z * Math.pow(10, 10);
			}
			content += "\n";

		}
		return content;
	}

	public static void recordPositions() {
		for (int i = 0; i < particles.size(); i++) {
			particles.get(i).recordPos();
		}
	}

	public static void moveParticlesGaus() {
		// moves/applies initial forces on the particles according to a gaussian
		// distribution
		for (int i = 0; i < particles.size(); i++) {
			particles.get(i).move(gau(gausMult), gau(gausMult), gau(gausMult));
		}

	}

	private static double getGausMult(double rms, boolean newWay) {
		if (constant) {
			return TMPgauMult;
		}

		//if (newWay) {
			//return ((double) ((2 * TMPgauMult * particleMass * frictionCoef * 1.3807 * Math.pow(10, -23) * targetTemp)
			//		/ Math.sqrt(3)));
		//}

		// calculates the gausMult that will yield highest probability of a translation
		// of (rms) units in 1d space
		// =====================Gaussian Multiplier Settings ===========================
		double GM = 1;
		double distanceSum = 0f;
		double avg = 0f;
		int sumIter = 100;
		int changeIter = 1000;
		double scaledRMS = scale * rms;
		// =============================================================================
		for (int j = 0; j < changeIter; j++) {
			distanceSum = 0f;
			avg = 0f;
			for (int i = 0; i < sumIter; i++) {
				distanceSum += Math.abs(gau(GM));
			}
			avg = distanceSum / sumIter;
			// adjust GM so avg gets closer to scaledrms
			if (avg > scaledRMS) {
				GM *= 1 - desentPerc;
			}
			if (avg < scaledRMS) {
				GM *= 1 + desentPerc;
			}
			// System.out.println(avg + ", " + GM);
		}

		return GM;
	}

	public static double gau(double _gausMult) {
		Random foo = new Random();
		return (double) (_gausMult * foo.nextGaussian());
	}

	// ======================= Functions ======================
	public static double arrayAvg(double[] arr) {
		double sum = 0f;
		for (int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}

		return (double) (sum / arr.length);

	}

	public static double alistAvg(ArrayList<Double> nums) {
		double sum = 0;
		for (int i = 0; i < nums.size(); i++) {
			sum += nums.get(i);
		}
		return (sum / nums.size());
	}

	public static extrema min(ArrayList<Double> nums) {
		extrema min = new extrema(0, nums.get(0));
		for (int i = 0; i < nums.size(); i++) {
			if (nums.get(i) < min.magnitude) {
				min = new extrema(i, nums.get(i));
			}
		}
		return min;
	}

	public static extrema max(ArrayList<Double> nums) {
		extrema max = new extrema(0, nums.get(0));
		for (int i = 0; i < nums.size(); i++) {
			if (nums.get(i) > max.magnitude) {
				max = new extrema(i, nums.get(i));
			}
		}
		return max;
	}

	public static void writeToFile(String path, String text) {
		// writes directly to file & will replace all previous text there
		try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(path))) {
			bufferedWriter.write(text);
		} catch (IOException e) {
			System.out.println("Error: IO Exception");
		}
	}

	public static Position sub(Position a, Position v) {
		return (new Position(a.x - v.x, a.y - v.y, a.z - v.z));
	}

	public static Position mult(Position pos, double v) {
		return (new Position(pos.x * v, pos.y * v, pos.z * v));
	}

	public static Position add(Position pos1, Position pos2) {
		return (new Position(pos1.x + pos2.x, pos1.y + pos2.y, pos1.z + pos2.z));
	}

	public static Position divide(Position pos, double num) {
		if (num == 0) {
			throw new IllegalArgumentException("Cannot Divide By 0");
		}
		double x = pos.x / num;
		double y = pos.y / num;
		double z = pos.z / num;

		return new Position(x, y, z);
	}

	public static double distance(Position pos1, Position pos2) {
		return (double) (Math.sqrt(((pos1.x - pos2.x) * (pos1.x - pos2.x)) + ((pos1.y - pos2.y) * (pos1.y - pos2.y))
				+ ((pos1.z - pos2.z) * (pos1.z - pos2.z))));
	}

	public static Position midpoint(Position pos1, Position pos2) {

		return (new Position((pos1.x + pos2.x) / 2, (pos1.y + pos2.y) / 2, (pos1.z + pos2.z) / 2));
	}

	public static Position vecSum(ArrayList<Vector> vecs) {
		double x = 0;
		double y = 0;
		double z = 0;
		for (int i = 0; i < vecs.size(); i++) {
			x += vecs.get(i).direction.x * vecs.get(i).magnitude;
			y += vecs.get(i).direction.y * vecs.get(i).magnitude;
			z += vecs.get(i).direction.z * vecs.get(i).magnitude;
		}
		return new Position(x, y, z);

	}

	public static double vecToMag(Position pos) {
		return ((Math.sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)));
	}

	public static double getRand() {
		Random rand = new Random();
		return ((rand.nextFloat() * 2) - 1);
	}

	public static double gradientFunction(double input) {
		return (Math.pow(2, input / substrateGradientCoeff));
	}

	public static double getMult(double input) {
		double temp = (((TMPMax - ((TMPMax + TMPLo) / 2)) / xBound) * input) + ((TMPMax + TMPLo) / 2);
		return (Math.sqrt(temp / (5.02 * (Math.pow(10, 18)))));
	}

}

class extrema {
	int timeStep;
	double magnitude;

	public extrema(int TtimeStep, double Tmagnitude) {
		timeStep = TtimeStep;
		magnitude = Tmagnitude;

	}
}

class Graph {

	public Graph(String title, int width, int height, ArrayList<Double> data) {
		JFrame frame = new JFrame(title);
		Canvas canvas = new Drawing(width, height, data);
		canvas.setSize(width, height);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(canvas);
		frame.pack();
		frame.setVisible(true);
	}
}