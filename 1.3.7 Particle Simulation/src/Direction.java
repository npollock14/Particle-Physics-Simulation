public class Direction {
	double x, y, z;

	public Direction(double x, double y, double z) {
		double mag = Math.sqrt((x * x) + (y * y) + (z * z));
		if (mag == 0) {
			this.x = 0;
			this.y = 0;
			this.z = 0;
		}else {
		this.x = x / mag;
		this.y = y / mag;
		this.z = z / mag;
		}
	}

	public double toMag() {
		return (Math.sqrt((x * x) + (y * y) + (z * z)));
	}
}