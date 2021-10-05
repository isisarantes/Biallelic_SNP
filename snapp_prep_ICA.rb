#!/usr/bin/env ruby

# Ísis Arantes, 2021-10-05
#
# An adaptation of Ruby script available at
# https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb
# to create files to run in the SNAPP (Beast2)
#
# This script prepares NEX format input files for the software SNAPP (http://beast2.org/snapp/),
# given a phylip format SNP matrix, a table linking species IDs and specimen IDs. 
#
# To learn about the available options of this script, run with '-h' (or without any options):
# ruby snapp_prep_ICA.rb -h
#
# For questions and bug fixes, email to isisarantes@gmail.com

# Load required libraries.
require 'optparse'

# Feedback.
puts ""
puts "snapp_prep_ICA.rb"
puts ""
puts "***************************************************"
puts ""

# Define default options.
options = {}
options[:phylip] = nil
options[:vcf] = nil
options[:table] = "example.spc.txt"
options[:max_snps] = nil
options[:transversions] = false
options[:transitions] = false
options[:nex] = "snapp.nex"
options[:out] = "snapp"
options[:no_annotation] = false

# Get the command line options.
ARGV << '-h' if ARGV.empty?
opt_parser = OptionParser.new do |opt|
	opt.banner = "Usage: ruby #{$0} [OPTIONS]"
	opt.separator  ""
	opt.separator  "Example:"
	opt.separator  "ruby #{$0} -p example.phy -t #{options[:table]} -x #{options[:nex]}"
	opt.separator  ""
	opt.separator  "Options:"
	opt.on("-p","--phylip FILENAME","File with SNP data in phylip format (default: none).") {|p| options[:phylip] = p}
	opt.on("-v","--vcf FILENAME","File with SNP data in vcf format (default: none).") {|v| options[:vcf] = v}
	opt.on("-t","--table FILENAME","File with table linking species and specimens (default: #{options[:table]}).") {|t| options[:table] = t}
	opt.on("-m","--max-snps NUMBER",Integer,"Maximum number of SNPs to be used (default: no maximum).") {|m| options[:max_snps] = m}
	opt.on("-r","--transversions","Use transversions only (default: #{options[:transversions]}).") {options[:transversions] = true}
	opt.on("-i","--transitions","Use transitions only (default: #{options[:transitions]}).") {options[:transitions] = true}
	opt.on("-x","--nex FILENAME","Output file in NEX format (default: #{options[:nex]}).") {|x| options[:nex] = x}
	opt.on("-h","--help","Print this help text.") {
		puts opt_parser
		exit(0)
	}
	opt.separator  ""
end
opt_parser.parse!

# Make sure that input is provided in either phylip or vcf format.
if options[:phylip] == nil and options[:vcf] == nil
	puts "ERROR: An input file must be provided, either in phylip format with option '-p' or in vcf format with option '-v'!"
	exit(1)
elsif options[:phylip] != nil and options[:vcf] != nil
	puts "ERROR: Only one of the two options '-p' and '-v' can be used!"
	exit(1)	
end

# Make sure that the -r and -i options are not used jointly.
if options[:transversions] and options[:transitions]
	puts "ERROR: Only one of the two options '-r' and '-i' can be used!"
	exit(1)
end

# Initiate a warn string and counts for excluded sites.
warn_string = ""
number_of_excluded_sites_missing = 0
number_of_excluded_sites_monomorphic = 0
number_of_excluded_sites_triallelic = 0
number_of_excluded_sites_tetraallelic = 0
number_of_excluded_sites_indel = 0
number_of_excluded_sites_transition = 0
number_of_excluded_sites_transversion = 0
number_of_sites_with_half_call = 0

# Initiate arrays for specimen ids and sequences.
specimen_ids = []
seqs = []

# Define various characters.
binary_chars = ["0","1","2"]
nucleotide_chars = ["A","C","G","T","R","Y","S","W","K","M"]
ambiguous_chars = ["R","Y","S","W","K","M"]
missing_chars = ["-","?","N"]

# Read the phylip file.
if options[:phylip] != nil
	phylip_file =  File.open(options[:phylip])
	phylip_lines = phylip_file.readlines
	phylip_lines[1..-1].each do |l|
		unless l.strip == ""
			specimen_ids << l.split[0]
			seqs << l.split[1].upcase
		end
	end
	# Recognize the sequence format (nucleotides or binary).
	unique_seq_chars = []
	seqs.each do |s|
		s.size.times do |pos|
			unless missing_chars.include?(s[pos])
				unique_seq_chars << s[pos] unless unique_seq_chars.include?(s[pos])
			end
		end
	end
	sequence_format_is_nucleotide = true
	sequence_format_is_binary = true
	unique_seq_chars.each do |c|
		sequence_format_is_binary = false unless binary_chars.include?(c)
		sequence_format_is_nucleotide = false unless nucleotide_chars.include?(c)
	end
	if sequence_format_is_binary == false and sequence_format_is_nucleotide == false
		puts "ERROR: Sequence format could not be recognized as either 'nucleotide' or 'binary'!"
		exit(1)
	end
end

# Alternatively, read the vcf file.
tmp_line_count_all = 0
if options[:vcf] != nil
	vcf_file =  File.open(options[:vcf])
	vcf_lines = vcf_file.readlines
	vcf_header_line = ""
	vcf_lines.each do |l|
		if l[0..1] != "##"
			if l[0] == "#" and vcf_header_line == ""
				vcf_header_line = l
				vcf_header_line_ary = vcf_header_line.split
				specimen_ids = vcf_header_line_ary[9..-1]
				specimen_ids.size.times {seqs << ""}
			elsif vcf_header_line != ""
				tmp_line_count_all += 1
				line_ary = l.split
				ref = line_ary[3]
				alt = line_ary[4]
				if ref.size == 1 and alt.size == 1
					format_str = line_ary[8]
					gt_index = format_str.split(":").index("GT")
					if gt_index == nil
						puts "ERROR: Expected 'GT' in FORMAT field but could not find it!"
						exit(1)
					end
					specimen_index = 0
					found_half_called_gt = false
					line_ary[9..-1].each do |rec|
						gt = rec.split(":")[gt_index]
						if gt.include?("/")
							gt1 = gt.split("/")[0]
							gt2 = gt.split("/")[1]
						elsif gt.include?("|")
							gt1 = gt.split("|")[0]
							gt2 = gt.split("|")[1]
						else
							puts "ERROR: Expected alleles to be separated by '/' or '|' but did not found such separators!"
							exit(1)
						end
						if ["0","1","."].include?(gt1) and ["0","1","."].include?(gt2)
							if gt1 == "0"
								base1 = ref
							elsif gt1 == "1"
								base1 = alt
							else
								base1 = "N"
							end
							if gt2 == "0"
								base2 = ref
							elsif gt2 == "1"
								base2 = alt
							else
								base2 = "N"
							end
							if base1 == "N" and base2 == "N"
								seqs[specimen_index] << "N"
							elsif base1 == "N" or base2 == "N"
								seqs[specimen_index] << "N"
								found_half_called_gt = true
							elsif [base1,base2].sort == ["A","A"]
								seqs[specimen_index] << "A"
							elsif [base1,base2].sort == ["A","C"]
								seqs[specimen_index] << "M"
							elsif [base1,base2].sort == ["A","G"]
								seqs[specimen_index] << "R"
							elsif [base1,base2].sort == ["A","T"]
								seqs[specimen_index] << "W"
							elsif [base1,base2].sort == ["C","C"]
								seqs[specimen_index] << "C"
							elsif [base1,base2].sort == ["C","G"]
								seqs[specimen_index] << "S"
							elsif [base1,base2].sort == ["C","T"]
								seqs[specimen_index] << "Y"
							elsif [base1,base2].sort == ["G","G"]
								seqs[specimen_index] << "G"
							elsif [base1,base2].sort == ["G","T"]
								seqs[specimen_index] << "K"
							elsif [base1,base2].sort == ["T","T"]
								seqs[specimen_index] << "T"
							else
								puts "ERROR: Unexpected genotype: #{base1} and #{base2}!"
								exit(1)
							end
						else
							puts "ERROR: Expected genotypes to be bi-allelic and contain only 0s and/or 1s or missing data marked with '.', but found #{gt1} and #{gt2}!"
							exit(1)
						end
						specimen_index += 1
					end
					number_of_sites_with_half_call += 1 if found_half_called_gt
				elsif ref.size > 1 and ref.include?(",") == false
					number_of_excluded_sites_indel += 1
				elsif alt.size > 1 and alt.include?(",") == false
					number_of_excluded_sites_indel += 1
				elsif ref.include?(",") or alt.include?(",")
					number_of_ref_and_alt_alleles = ref.count(",") + alt.count(",") + 2
					if number_of_ref_and_alt_alleles == 3
						number_of_excluded_sites_triallelic += 1
					elsif number_of_ref_and_alt_alleles == 4
						number_of_excluded_sites_tetraallelic += 1
					else
						puts "ERROR: Unexpected combination of REF and ALT alleles (REF: #{ref}; ALT: #{alt})!"
						exit(1)
					end
				end
			else
				puts "ERROR: Expected a vcf header line beginning with '#CHROM' but could not find it!"
				exit(1)
			end
		end
	end

	# Make sure that all sequences have a positive and equal length.
	seqs[1..-1].each do |s|
		if s.size != seqs[0].size
			puts "ERROR: Sequences have different lengths!"
			exit(1)
		end
	end

	# Specify that the sequence format is nucleotide.
	sequence_format_is_binary = false
end

# If necessary, translate the sequences into SNAPP's "0", "1", "2" code, where "1" is heterozygous.
binary_seqs = []
seqs.size.times{binary_seqs << ""}
if sequence_format_is_binary
	seqs[0].size.times do |pos|
		alleles_at_this_pos = []
		seqs.each do |s|
			alleles_at_this_pos << s[pos] if binary_chars.include?(s[pos])
		end
		uniq_alleles_at_this_pos = alleles_at_this_pos.uniq
		if uniq_alleles_at_this_pos.size == 2 or uniq_alleles_at_this_pos.size == 3
			seqs.size.times do |x|
				binary_seqs[x] << seqs[x][pos]
			end
		elsif uniq_alleles_at_this_pos.size == 0
			number_of_excluded_sites_missing += 1
		elsif uniq_alleles_at_this_pos.size == 1
			number_of_excluded_sites_monomorphic += 1
		end
	end

	# Check if the total number of '0' and '2' in the data set are similar.
	total_number_of_0s = 0
	total_number_of_2s = 0
	binary_seqs.each do |b|
		total_number_of_0s += b.count("0")
		total_number_of_2s += b.count("2")
	end
	goal = 0.5
	tolerance = 0.01
	if (total_number_of_0s/(total_number_of_0s+total_number_of_2s).to_f - goal).abs > tolerance
		warn_string << "WARNING: The number of '0' and '2' in the data set is expected to be similar, however,\n"
		warn_string << "    they differ by more than #{tolerance*100.round} percent.\n"
		warn_string << "\n"
	end
else
	seqs[0].size.times do |pos|
		# Collect all bases at this position.
		bases_at_this_pos = []
		seqs.each do |s|
			if s[pos] == "A"
				bases_at_this_pos << "A"
				bases_at_this_pos << "A"
			elsif s[pos] == "C"
				bases_at_this_pos << "C"
				bases_at_this_pos << "C"
			elsif s[pos] == "G"
				bases_at_this_pos << "G"
				bases_at_this_pos << "G"
			elsif s[pos] == "T"
				bases_at_this_pos << "T"
				bases_at_this_pos << "T"
			elsif s[pos] == "R"
				bases_at_this_pos << "A"
				bases_at_this_pos << "G"
			elsif s[pos] == "Y"
				bases_at_this_pos << "C"
				bases_at_this_pos << "T"
			elsif s[pos] == "S"
				bases_at_this_pos << "G"
				bases_at_this_pos << "C"
			elsif s[pos] == "W"
				bases_at_this_pos << "A"
				bases_at_this_pos << "T"
			elsif s[pos] == "K"
				bases_at_this_pos << "G"
				bases_at_this_pos << "T"
			elsif s[pos] == "M"
				bases_at_this_pos << "A"
				bases_at_this_pos << "C"
			else
				unless ["N","?","-"].include?(s[pos])
					puts "ERROR: Found unexpected base at position #{pos+1}: #{s[pos]}!"
					exit(1)
				end
			end
		end
		uniq_bases_at_this_pos = bases_at_this_pos.uniq
		# Issue a warning if non-bi-allelic sites are excluded.
		if uniq_bases_at_this_pos.size == 2
			# Determine if this is a transition or transversion site.
			transversion_site = false
			if uniq_bases_at_this_pos.sort == ["A","C"]
				transversion_site = true
			elsif uniq_bases_at_this_pos.sort == ["A","G"]
				transversion_site = false
			elsif uniq_bases_at_this_pos.sort == ["A","T"]
				transversion_site = true
			elsif uniq_bases_at_this_pos.sort == ["C","G"]
				transversion_site = true
			elsif uniq_bases_at_this_pos.sort == ["C","T"]
				transversion_site = false
			elsif uniq_bases_at_this_pos.sort == ["G","T"]
				transversion_site = true
			else
				puts "ERROR: Unexpected combination of unique bases at position #{pos+1}: #{uniq_bases_at_this_pos[0]}, #{uniq_bases_at_this_pos[1]}"
				exit(1)
			end
			transition_site = true
			transition_site = false if transversion_site == true
			if options[:transversions] == true and transversion_site == false # If the site is a transition and only transversions are allowed.
				number_of_excluded_sites_transition += 1
			elsif options[:transitions] == true and transition_site == false # If the site is a transversion and only transitions are allowed.
				number_of_excluded_sites_transversion += 1
			else
				# Randomly define what's "0" and "2".
				uniq_bases_at_this_pos.shuffle!
				seqs.size.times do |x|
					if seqs[x][pos] == uniq_bases_at_this_pos[0]
						binary_seqs[x] << "0"
					elsif seqs[x][pos] == uniq_bases_at_this_pos[1]
						binary_seqs[x] << "2"
					elsif missing_chars.include?(seqs[x][pos])
						binary_seqs[x] << "-"
					elsif ambiguous_chars.include?(seqs[x][pos])
						binary_seqs[x] << "1"
					else
						puts "ERROR: Found unexpected base at position #{pos+1}: #{seqs[x][pos]}!"
						exit(1)
					end
				end
			end
		elsif uniq_bases_at_this_pos.size == 0
			number_of_excluded_sites_missing += 1
		elsif uniq_bases_at_this_pos.size == 1
			number_of_excluded_sites_monomorphic += 1
		elsif uniq_bases_at_this_pos.size == 3
			number_of_excluded_sites_triallelic += 1
		elsif uniq_bases_at_this_pos.size == 4
			number_of_excluded_sites_tetraallelic += 1
		else
			puts "ERROR: Found unexpected number of alleles at position #{pos+1}!"
			exit(1)
		end
	end
end

# Read the file with a table linking species and specimens.
table_file = File.open(options[:table])
table_lines = table_file.readlines
table_species = []
table_specimens = []
table_lines.each do |l|
	line_ary = l.split
	header_line = false
	header_line = true if line_ary[0].downcase == "species" and line_ary[1].downcase == "specimen"
	header_line = true if line_ary[0].downcase == "species" and line_ary[1].downcase == "specimens"
	header_line = true if line_ary[0].downcase == "species" and line_ary[1].downcase == "sample"
	header_line = true if line_ary[0].downcase == "species" and line_ary[1].downcase == "samples"
	unless header_line
		table_species << line_ary[0]
		table_specimens << line_ary[1]
	end
end

# Make sure that the arrays table_specimens and specimen_ids are identical when sorted.
unless table_specimens.sort == specimen_ids.sort
	puts "ERROR: The specimens listed in file #{options[:table]} and those included in the input file are not identical!"
	exit(1)
end

# Remove sites at which one or more species have only missing data; these could not be used by SNAPP anyway.
binary_seqs_for_snapp = []
binary_seqs.size.times {binary_seqs_for_snapp << ""}
binary_seqs[0].size.times do |pos|
	one_or_more_species_have_only_missing_data_at_this_pos = false
	table_species.uniq.each do |spc|
		specimens_for_this_species = []
		table_specimens.size.times do |x|
			specimens_for_this_species << table_specimens[x] if table_species[x] == spc
		end
		alleles_for_this_species_at_this_pos = []
		specimen_ids.size.times do |x|
			if specimens_for_this_species.include?(specimen_ids[x])
				alleles_for_this_species_at_this_pos << binary_seqs[x][pos]
			end
		end
		if alleles_for_this_species_at_this_pos.uniq == ["-"]
			one_or_more_species_have_only_missing_data_at_this_pos =  true
		end
	end
	# Set all alleles at this position to nil if one species had only missing data.
	if one_or_more_species_have_only_missing_data_at_this_pos
		number_of_excluded_sites_missing += 1
	else
		binary_seqs.size.times do |x|
			binary_seqs_for_snapp[x] << binary_seqs[x][pos]
		end
	end
end
binary_seqs = binary_seqs_for_snapp

# If a maximum number of SNPs has been set, reduce the data set to this number.
number_of_sites_before_excluding_due_to_max = binary_seqs[0].size
number_of_excluded_sites_due_to_max = 0
if options[:max_snps] != nil
	if options[:max_snps] < binary_seqs[0].size
		seq_indices = []
		binary_seqs[0].size.times {|x| seq_indices << x}
		selected_seq_indices = seq_indices.sample(options[:max_snps]).sort
		binary_seqs_red = []
		binary_seqs.each do |s|
			binary_seq_red = ""
			selected_seq_indices.each do |i|
				binary_seq_red << s[i]
			end
			binary_seqs_red << binary_seq_red
		end
		binary_seqs = binary_seqs_red
		number_of_excluded_sites_due_to_max = number_of_sites_before_excluding_due_to_max - options[:max_snps]
	else
		warn_string << "WARNING: The maximum number of SNPs has been set to #{options[:max_snps]}, which is greater\n"
		warn_string << "    than the number of bi-allelic SNPs with sufficient information for SNAPP.\n"
	end
end

# Compose the warn string if necessary.
if number_of_sites_with_half_call > 0
	warn_string << "WARNING: Found #{number_of_sites_with_half_call} site"
	warn_string << "s" if number_of_sites_with_half_call > 1
	warn_string << " with genotypes that were half missing. These genotypes were ignored.\n"
end
if number_of_excluded_sites_missing > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_missing} site"
	warn_string << "s" if number_of_excluded_sites_missing > 1
	warn_string << " with only missing data in one or more species.\n"
end
if number_of_excluded_sites_monomorphic > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_monomorphic} monomorphic site"
	warn_string << "s" if number_of_excluded_sites_monomorphic > 1
	warn_string << ".\n"
end
if number_of_excluded_sites_transition > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_transition} transition site"
	warn_string << "s" if number_of_excluded_sites_transition > 1
	warn_string << ".\n"
end
if number_of_excluded_sites_transversion > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_transversion} transversion site"
	warn_string << "s" if number_of_excluded_sites_transversion > 1
	warn_string << ".\n"
end
if number_of_excluded_sites_triallelic > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_triallelic} tri-allelic site"
	warn_string << "s" if number_of_excluded_sites_triallelic > 1
	warn_string << ".\n"
end
if number_of_excluded_sites_tetraallelic > 0
	warn_string << "WARNING: Excluded #{number_of_excluded_sites_tetraallelic} tetra-allelic site"
	warn_string << "s" if number_of_excluded_sites_tetraallelic > 1
	warn_string << ".\n"
end

# If there were any warning, print them.
unless warn_string == ""
	warn_string << "\n"
	puts warn_string
end

# Print the info string.
if options[:max_snps] != nil
	info_string = "INFO: Removed #{number_of_excluded_sites_due_to_max} bi-allelic sites due to specified maximum number of #{options[:max_snps]} sites.\n"
	info_string << "\n"
	puts info_string
else
	if options[:transversions]
		info_string = "INFO: Retained #{binary_seqs[0].size} bi-allelic transversion sites.\n"
	elsif options[:transitions]
		info_string = "INFO: Retained #{binary_seqs[0].size} bi-allelic transition sites.\n"
	else
		info_string = "INFO: Retained #{binary_seqs[0].size} bi-allelic sites.\n"
	end
	info_string << "\n"
	puts info_string
end


# Prepare SNAPP input string.
snapp_string = ""
snapp_string << "#NEXUS\n"
snapp_string << "\n"
unless options[:no_annotation]
	snapp_string << "\n"
	if sequence_format_is_binary
		snapp_string << "[The SNP data matrix, taken from file #{options[:phylip]}.]\n"
	else
		snapp_string << "[The SNP data matrix, converted to binary format from file #{options[:phylip]}.]\n"
	end
	snapp_string << "\n"
end
snapp_string << "Begin data;\n"
snapp_string << "	Dimensions ntax=#{table_specimens.length} nchar=#{binary_seqs[0].size};\n"
snapp_string <<	"	Format datatype=integerdata symbols='012' gap=-;\n"
snapp_string << "	Matrix\n"
specimen_ids.size.times do |x|
	snapp_string << "#{specimen_ids[x]}\_#{table_species[x]}	#{binary_seqs[x]}\n"
end
snapp_string << "	;\n"
snapp_string << "End;\n"

# Write the SNAPP input file.
snapp_file = File.open(options[:nex],"w")
snapp_file.write(snapp_string)
puts "Wrote SNAPP input in NEX format to file #{options[:nex]}.\n\n"
