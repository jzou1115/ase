package genome;

import java.util.Iterator;

import snp.SNP;

public class GenomicRegion {
	GenomicCoordinate start;
	GenomicCoordinate end;

	public GenomicRegion(GenomicCoordinate s, GenomicCoordinate e){
		start=s;
		end=e;
	}
	
	public long size(){
		return start.distance(end) + 1;
	}
	
	public int getChromosome(){
		return start.getChromosome();
	}
	
	public GenomicCoordinate getStart(){
		return start;
	}
	
	public GenomicCoordinate getEnd(){
		return end;
	}
	
	public GenomicRegion increment(int numBases){
		return new GenomicRegion(start.increment(numBases), end.increment(numBases));
	}
	
	public GenomicRegion decrement(int numBases){
		return new GenomicRegion(start.decrement(numBases), end.decrement(numBases));
	}
	
	//make sure this has genomic relevance
	public GenomicRegion expand(int numBases){
		return new GenomicRegion(start.decrement(numBases), end.increment(numBases));		
	}
	
	public int toIndex(GenomicCoordinate c){
		if(!this.contains(c)){
			throw new RuntimeException("Region "+this+" does not contain coordinate "+c+".");
		}
		return (int) start.distance(c);
	}

	public GenomicCoordinate toCoordinate(int index){
		if(index >= this.size()){
			throw new RuntimeException("Index "+index+" out of bounds on region "+this);
		}
		return start.increment(index);
	}
	
	public long distance(GenomicRegion other){
		if(this.overlaps(other)){
			return -this.getOverlap(other);
		}
		if(start.compareTo(other.getStart()) < 0){
			return end.distance(other.getStart());
		}
		return other.getEnd().distance(start);
	}
	
	public boolean contains(GenomicCoordinate coordinate){
		return start.compareTo(coordinate) <= 0 && end.compareTo(coordinate) >= 0;
	}
	
	
	public boolean contains(GenomicRegion other){
		return this.contains(other.start) && this.contains(other.end);
	}
	
	public boolean overlaps(GenomicRegion other){
		return other.contains(this) || this.contains(other.start) || this.contains(other.end);
	}
	
	public long getOverlap(GenomicRegion other){
		GenomicRegion overlap = this.intersection(other);
		if(overlap == null) return 0;
		return overlap.size();
	}
	
	public GenomicRegion intersection(GenomicRegion other){
		if(!this.overlaps(other)){
			return null;
		}
		if(this.contains(other)){
			return other;
		}
		if(other.contains(this)){
			return this;
		}
		if(this.contains(other.start)){
			return new GenomicRegion(other.start, end);
		}
		return new GenomicRegion(start, other.end);
	}
	
	
	
	/**
	 * Splits this genomic region into two regions around the given coordinate. The left region
	 * is [start, coordinate) and the right region is [coordinate, end]. This will error
	 * if the coordinate is not within this region or the coordinate is defined using a different
	 * genome from this region.
	 * <p>
	 * Note that because the coordinate is included in the right side region, if the coordinate is
	 * the start of this region, then the right side region will be the full region and the left
	 * side region will contain no bases. Because a region cannot contain zero bases, the full region
	 * will be returned as BOTH the left and right side regions of the array.
	 * @param coordinate - around which this region should be split
	 * @return an array of {left region, right region}
	 */
	public GenomicRegion[] split(GenomicCoordinate coordinate){
		if(!this.contains(coordinate)){
			throw new RuntimeException("The region "+this+" cannot be split around the coordinate "+coordinate);
		}
		if(start.equals(coordinate)){
			return new GenomicRegion[]{this, this};
		}
		GenomicRegion left = new GenomicRegion(start, coordinate.decrement(1));
		GenomicRegion right = new GenomicRegion(coordinate, end);
		return new GenomicRegion[]{left, right};
	}
	
	@Override
	public String toString(){
		return start.getChromosome()+"\t"+start.getCoord()+"\t"+end.getCoord();
	}
	
	@Override
	public boolean equals(Object o){
		if(o == null) return false;
		if(this == o) return true;
		if(o instanceof GenomicRegion){
			GenomicRegion other = (GenomicRegion) o;
			return start.equals(other.start) && end.equals(other.end);
		}
		return false;
	}
	
	public int compareTo(GenomicRegion other) {
		return this.getStart().compareTo(other.getStart());
	}
	
	@Override
	public int hashCode(){
		return (int) (7*start.getCoord()+13*end.getCoord());
	}


	
}
