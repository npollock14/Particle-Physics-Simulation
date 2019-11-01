public class Position {
	public double x;
	public double y;
	public double z;

	public Position(double Tx, double Ty, double Tz) {
		x = Tx;
		y = Ty;
		z = Tz;
	}

	public Direction toDir() {
		return (new Direction(this.x, this.y, this.z));
	}

	public double toMag() {
		return (Math.sqrt((x * x) + (y * y) + (z * z)));
	}
}