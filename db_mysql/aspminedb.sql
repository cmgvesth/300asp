-- phpMyAdmin SQL Dump
-- version 4.0.10deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost:3306
-- Generation Time: Feb 18, 2015 at 01:52 PM
-- Server version: 5.5.41-0ubuntu0.14.04.1
-- PHP Version: 5.5.9-1ubuntu4.5

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `aspminedb`
--
CREATE DATABASE IF NOT EXISTS `aspminedb` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `aspminedb`;

-- --------------------------------------------------------

--
-- Table structure for table `antiblast`
--

DROP TABLE IF EXISTS `antiblast`;
CREATE TABLE IF NOT EXISTS `antiblast` (
  `org_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `clust_id` varchar(74) NOT NULL DEFAULT '',
  `name` varchar(100) NOT NULL,
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL,
  `max_pident` decimal(10,2) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `antismash`
--

DROP TABLE IF EXISTS `antismash`;
CREATE TABLE IF NOT EXISTS `antismash` (
  `org_id` int(11) NOT NULL,
  `filename` text NOT NULL,
  `clust_size` int(11) NOT NULL,
  `clust_backbone` varchar(50) NOT NULL,
  `clust_origin` varchar(50) NOT NULL,
  `clust_start` int(11) NOT NULL,
  `clust_end` int(11) NOT NULL,
  `clust_type` varchar(50) NOT NULL,
  `sm_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `sm_short` varchar(100) NOT NULL,
  `sm_desc` text NOT NULL,
  `cluster` text NOT NULL,
  PRIMARY KEY (`clust_backbone`,`sm_protein_id`),
  KEY `sm_protein_id` (`sm_protein_id`),
  KEY `clust_backbone` (`clust_backbone`),
  KEY `clust_backbone_sm_protein_id` (`clust_backbone`,`sm_protein_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `antismash2blast`
--

DROP TABLE IF EXISTS `antismash2blast`;
CREATE TABLE IF NOT EXISTS `antismash2blast` (
  `org_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `clust_id` varchar(74) NOT NULL DEFAULT '',
  `name` varchar(100) NOT NULL,
  `q_org` varchar(100) NOT NULL,
  `q_seqkey` varchar(100) NOT NULL,
  `q_tail` varchar(100) NOT NULL,
  `h_org` varchar(100) NOT NULL,
  `h_seqkey` varchar(100) NOT NULL,
  `h_tail` varchar(100) NOT NULL,
  `pident` decimal(10,2) DEFAULT NULL,
  `q_len` int(11) NOT NULL,
  `q_start` int(11) NOT NULL,
  `q_end` int(11) NOT NULL,
  `h_len` int(11) NOT NULL,
  `h_start` int(11) NOT NULL,
  `h_end` int(11) NOT NULL,
  `bitscore` int(11) NOT NULL,
  `evalue` varchar(10) NOT NULL,
  `h_seqid` varchar(100) NOT NULL,
  `q_seqid` varchar(100) NOT NULL,
  `filename` varchar(100) NOT NULL,
  `loadstamp` varchar(100) NOT NULL,
  `q_cov` decimal(10,2) DEFAULT NULL,
  `h_cov` decimal(10,2) DEFAULT NULL,
  KEY `sm_protein_id` (`sm_protein_id`,`name`),
  KEY `name` (`name`),
  KEY `sm_protein_id_2` (`sm_protein_id`),
  KEY `seqids` (`h_seqid`,`q_seqid`),
  KEY `qseq_horg_pident` (`q_seqid`,`h_org`,`pident`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `antismash2organism`
--

DROP TABLE IF EXISTS `antismash2organism`;
CREATE TABLE IF NOT EXISTS `antismash2organism` (
  `org_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `clust_id` varchar(74) NOT NULL DEFAULT '',
  `name` varchar(100) NOT NULL,
  KEY `sm_protein_id` (`sm_protein_id`,`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `antitmp`
--

DROP TABLE IF EXISTS `antitmp`;
CREATE TABLE IF NOT EXISTS `antitmp` (
  `org_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `clust_id` varchar(86) NOT NULL DEFAULT '',
  `name` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `assembly`
--

DROP TABLE IF EXISTS `assembly`;
CREATE TABLE IF NOT EXISTS `assembly` (
  `assembly_orgkey` varchar(100) NOT NULL,
  `assembly_seqkey` varchar(100) NOT NULL,
  `assembly_tailkey` varchar(100) NOT NULL,
  `org_id` int(11) NOT NULL,
  `assembly_seqname` varchar(100) NOT NULL,
  `assembly_orgseq` varchar(100) DEFAULT NULL,
  `assembly_seq` longblob,
  PRIMARY KEY (`assembly_seqname`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `blast`
--

DROP TABLE IF EXISTS `blast`;
CREATE TABLE IF NOT EXISTS `blast` (
  `q_org` varchar(100) NOT NULL,
  `q_seqkey` varchar(100) NOT NULL,
  `q_tail` varchar(100) NOT NULL,
  `h_org` varchar(100) NOT NULL,
  `h_seqkey` varchar(100) NOT NULL,
  `h_tail` varchar(100) NOT NULL,
  `pident` decimal(10,2) DEFAULT NULL,
  `q_len` int(11) NOT NULL,
  `q_start` int(11) NOT NULL,
  `q_end` int(11) NOT NULL,
  `h_len` int(11) NOT NULL,
  `h_start` int(11) NOT NULL,
  `h_end` int(11) NOT NULL,
  `bitscore` int(11) NOT NULL,
  `evalue` varchar(10) NOT NULL,
  `h_seqid` varchar(100) NOT NULL,
  `q_seqid` varchar(100) NOT NULL,
  `filename` varchar(100) NOT NULL,
  `loadstamp` varchar(100) NOT NULL,
  `q_cov` decimal(10,2) DEFAULT NULL,
  `h_cov` decimal(10,2) DEFAULT NULL,
  PRIMARY KEY (`bitscore`,`q_org`,`q_seqkey`,`h_org`,`h_seqkey`,`q_start`,`h_start`),
  KEY `q_org` (`q_org`),
  KEY `q_seqkey` (`q_seqkey`),
  KEY `h_org` (`h_org`),
  KEY `h_seqkey` (`h_seqkey`),
  KEY `h_seqid` (`h_seqid`),
  KEY `q_seqid` (`q_seqid`),
  KEY `filename` (`filename`),
  KEY `q_seqkey_2` (`q_seqkey`,`q_org`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `cds`
--

DROP TABLE IF EXISTS `cds`;
CREATE TABLE IF NOT EXISTS `cds` (
  `org_id` int(11) NOT NULL,
  `cds_orgseq` varchar(100) DEFAULT NULL,
  `cds_orgkey` varchar(100) NOT NULL,
  `cds_seqkey` varchar(100) NOT NULL,
  `cds_tailkey` varchar(100) NOT NULL,
  `cds_seqname` varchar(100) NOT NULL,
  `cds_seq` longblob,
  PRIMARY KEY (`cds_seqname`,`org_id`),
  KEY `org_id` (`org_id`),
  KEY `icds_tailkey` (`cds_tailkey`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `environment`
--

DROP TABLE IF EXISTS `environment`;
CREATE TABLE IF NOT EXISTS `environment` (
  `env_id` int(11) NOT NULL AUTO_INCREMENT,
  `env_desc` varchar(100) NOT NULL,
  `env_class` varchar(100) NOT NULL,
  PRIMARY KEY (`env_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `gff`
--

DROP TABLE IF EXISTS `gff`;
CREATE TABLE IF NOT EXISTS `gff` (
  `gff_seqorigin` varchar(100) NOT NULL,
  `org_id` int(11) NOT NULL DEFAULT '0',
  `gff_name` varchar(100) NOT NULL,
  `gff_type` varchar(100) NOT NULL,
  `gff_protein_id` varchar(100) DEFAULT NULL,
  `gff_trans_id` varchar(100) DEFAULT NULL,
  `gff_attributes` varchar(500) NOT NULL DEFAULT '',
  `gff_start` int(11) NOT NULL DEFAULT '0',
  `gff_end` int(11) DEFAULT NULL,
  `gff_length` int(11) DEFAULT NULL,
  `gff_score` varchar(10) DEFAULT NULL,
  `gff_strand` varchar(10) DEFAULT NULL,
  `gff_phase` varchar(10) DEFAULT NULL,
  `gff_exonnr` varchar(500) NOT NULL,
  PRIMARY KEY (`org_id`,`gff_attributes`,`gff_start`),
  KEY `igff_seq_id` (`gff_protein_id`),
  KEY `igff_name` (`gff_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `go`
--

DROP TABLE IF EXISTS `go`;
CREATE TABLE IF NOT EXISTS `go` (
  `go_term_id` varchar(50) NOT NULL DEFAULT '',
  `go_name` varchar(500) DEFAULT NULL,
  `go_termtype` varchar(50) DEFAULT NULL,
  `go_acc` varchar(50) NOT NULL DEFAULT '',
  PRIMARY KEY (`go_term_id`,`go_acc`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ipr`
--

DROP TABLE IF EXISTS `ipr`;
CREATE TABLE IF NOT EXISTS `ipr` (
  `ipr_id` varchar(50) NOT NULL DEFAULT '',
  `ipr_desc` text,
  `ipr_domaindb` varchar(50) DEFAULT NULL,
  `ipr_domain_id` varchar(100) NOT NULL DEFAULT '',
  `ipr_domaindesc` text,
  PRIMARY KEY (`ipr_id`,`ipr_domain_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `kegg`
--

DROP TABLE IF EXISTS `kegg`;
CREATE TABLE IF NOT EXISTS `kegg` (
  `kegg_id` int(11) NOT NULL AUTO_INCREMENT,
  `kegg_ecNum` varchar(50) NOT NULL,
  `kegg_definition` varchar(500) DEFAULT NULL,
  `kegg_catalyticActivity` varchar(500) DEFAULT NULL,
  `kegg_cofactors` varchar(500) DEFAULT NULL,
  `kegg_associatedDiseases` varchar(500) DEFAULT NULL,
  `kegg_pathway` varchar(500) NOT NULL,
  `kegg_pathway_class` varchar(500) DEFAULT NULL,
  `kegg_pathway_type` varchar(500) DEFAULT NULL,
  PRIMARY KEY (`kegg_id`),
  UNIQUE KEY `kegg_ecNum_kegg_pathway` (`kegg_ecNum`,`kegg_pathway`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=406029 ;

-- --------------------------------------------------------

--
-- Table structure for table `kog`
--

DROP TABLE IF EXISTS `kog`;
CREATE TABLE IF NOT EXISTS `kog` (
  `kog_id` varchar(50) NOT NULL DEFAULT '',
  `kog_defline` varchar(500) DEFAULT NULL,
  `kog_Class` varchar(500) DEFAULT NULL,
  `kog_Group` varchar(500) DEFAULT NULL,
  PRIMARY KEY (`kog_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lengths`
--

DROP TABLE IF EXISTS `lengths`;
CREATE TABLE IF NOT EXISTS `lengths` (
  `org_id` int(11) NOT NULL,
  `prot_seqname` varchar(100) NOT NULL,
  `len` int(10) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_extended`
--

DROP TABLE IF EXISTS `map_niger_extended`;
CREATE ALGORITHM=UNDEFINED DEFINER=`asp`@`localhost` SQL SECURITY DEFINER VIEW `map_niger_extended` AS (select `aspminedb`.`map_niger`.`query_id` AS `query_id`,(case when (`aspminedb`.`map_niger`.`hit_org` = 'Aniger1') then `aspminedb`.`map_niger`.`hit_id` end) AS `Aniger1`,(case when (`aspminedb`.`map_niger`.`hit_org` = 'Aniger3') then `aspminedb`.`map_niger`.`hit_id` end) AS `Aniger3`,(case when (`aspminedb`.`map_niger`.`hit_org` = 'Aniger7') then `aspminedb`.`map_niger`.`hit_id` end) AS `Aniger7` from `map_niger`);

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_extended_Pivot`
--

DROP TABLE IF EXISTS `map_niger_extended_Pivot`;
CREATE ALGORITHM=UNDEFINED DEFINER=`asp`@`localhost` SQL SECURITY DEFINER VIEW `map_niger_extended_Pivot` AS (select `map_niger_extended`.`query_id` AS `query_id`,`map_niger_extended`.`Aniger1` AS `Aniger1`,`map_niger_extended`.`Aniger3` AS `Aniger3`,`map_niger_extended`.`Aniger7` AS `Aniger7` from `map_niger_extended` group by `map_niger_extended`.`query_id`);

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_prot`
--

DROP TABLE IF EXISTS `map_niger_prot`;
CREATE TABLE IF NOT EXISTS `map_niger_prot` (
  `filename` varchar(50) DEFAULT NULL,
  `hit_id` varchar(100) DEFAULT NULL,
  `hit_org` varchar(100) DEFAULT NULL,
  `query_id` varchar(100) DEFAULT NULL,
  `query_org` varchar(100) DEFAULT NULL,
  `query_range` varchar(11) DEFAULT NULL,
  `hit_range` varchar(11) DEFAULT NULL,
  `aln_span` int(11) DEFAULT NULL,
  `alphabet` varchar(11) DEFAULT NULL,
  `bitscore` decimal(10,2) DEFAULT NULL,
  `bitscore_raw` decimal(10,2) DEFAULT NULL,
  `evalue` varchar(11) DEFAULT NULL,
  `gap_num` int(11) DEFAULT NULL,
  `hit_end` int(11) DEFAULT NULL,
  `hit_frame` varchar(11) DEFAULT NULL,
  `hit_span` int(11) DEFAULT NULL,
  `hit_start` int(11) DEFAULT NULL,
  `hit_strand` int(11) DEFAULT NULL,
  `ident_num` int(11) DEFAULT NULL,
  `is_fragmented` varchar(11) DEFAULT NULL,
  `pos_num` int(11) DEFAULT NULL,
  `query_end` int(11) DEFAULT NULL,
  `query_frame` varchar(11) DEFAULT NULL,
  `query_span` int(11) DEFAULT NULL,
  `query_start` int(11) DEFAULT NULL,
  `query_strand` int(11) DEFAULT NULL,
  `pident` decimal(10,2) DEFAULT NULL,
  UNIQUE KEY `filename` (`filename`,`hit_id`,`query_id`,`bitscore_raw`,`pident`),
  KEY `hit_id` (`hit_id`),
  KEY `query_id` (`query_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_prot_complete`
--

DROP TABLE IF EXISTS `map_niger_prot_complete`;
CREATE TABLE IF NOT EXISTS `map_niger_prot_complete` (
  `query_id` varchar(100) DEFAULT NULL,
  `query_org` varchar(100) DEFAULT NULL,
  `A1` varchar(100) DEFAULT NULL,
  `A1cov` decimal(14,0) DEFAULT NULL,
  `A1p` decimal(9,0) DEFAULT NULL,
  `A3` varchar(100) DEFAULT NULL,
  `A3cov` decimal(14,0) DEFAULT NULL,
  `A3p` decimal(9,0) DEFAULT NULL,
  `A7` varchar(100) DEFAULT NULL,
  `A7cov` decimal(14,0) DEFAULT NULL,
  `A7p` decimal(9,0) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_trans`
--

DROP TABLE IF EXISTS `map_niger_trans`;
CREATE TABLE IF NOT EXISTS `map_niger_trans` (
  `filename` varchar(50) DEFAULT NULL,
  `hit_id` varchar(100) DEFAULT NULL,
  `hit_org` varchar(100) DEFAULT NULL,
  `query_id` varchar(100) DEFAULT NULL,
  `query_org` varchar(100) DEFAULT NULL,
  `query_range` varchar(11) DEFAULT NULL,
  `hit_range` varchar(11) DEFAULT NULL,
  `aln_span` int(11) DEFAULT NULL,
  `alphabet` varchar(11) DEFAULT NULL,
  `bitscore` decimal(10,2) DEFAULT NULL,
  `bitscore_raw` decimal(10,2) DEFAULT NULL,
  `evalue` varchar(11) DEFAULT NULL,
  `gap_num` int(11) DEFAULT NULL,
  `hit_end` int(11) DEFAULT NULL,
  `hit_frame` varchar(11) DEFAULT NULL,
  `hit_span` int(11) DEFAULT NULL,
  `hit_start` int(11) DEFAULT NULL,
  `hit_strand` int(11) DEFAULT NULL,
  `ident_num` int(11) DEFAULT NULL,
  `is_fragmented` varchar(11) DEFAULT NULL,
  `pos_num` int(11) DEFAULT NULL,
  `query_end` int(11) DEFAULT NULL,
  `query_frame` varchar(11) DEFAULT NULL,
  `query_span` int(11) DEFAULT NULL,
  `query_start` int(11) DEFAULT NULL,
  `query_strand` int(11) DEFAULT NULL,
  `pident` decimal(10,2) DEFAULT NULL,
  UNIQUE KEY `filename` (`filename`,`hit_id`,`query_id`,`bitscore_raw`,`pident`),
  KEY `hit_id` (`hit_id`),
  KEY `query_id` (`query_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `map_niger_trans_complete`
--

DROP TABLE IF EXISTS `map_niger_trans_complete`;
CREATE TABLE IF NOT EXISTS `map_niger_trans_complete` (
  `query_id` varchar(100) DEFAULT NULL,
  `query_org` varchar(100) DEFAULT NULL,
  `A1` varchar(100) DEFAULT NULL,
  `A1cov` decimal(14,0) DEFAULT NULL,
  `A1p` decimal(9,0) DEFAULT NULL,
  `A3` varchar(100) DEFAULT NULL,
  `A3cov` decimal(14,0) DEFAULT NULL,
  `A3p` decimal(9,0) DEFAULT NULL,
  `A7` varchar(100) DEFAULT NULL,
  `A7cov` decimal(14,0) DEFAULT NULL,
  `A7p` decimal(9,0) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `organism`
--

DROP TABLE IF EXISTS `organism`;
CREATE TABLE IF NOT EXISTS `organism` (
  `org_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `real_name` text,
  `section` varchar(100) NOT NULL,
  `status` varchar(50) NOT NULL,
  PRIMARY KEY (`org_id`),
  KEY `iname` (`name`),
  KEY `lable` (`section`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=52 ;

-- --------------------------------------------------------

--
-- Table structure for table `org_has_environment`
--

DROP TABLE IF EXISTS `org_has_environment`;
CREATE TABLE IF NOT EXISTS `org_has_environment` (
  `env_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  PRIMARY KEY (`env_id`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `org_is_species`
--

DROP TABLE IF EXISTS `org_is_species`;
CREATE TABLE IF NOT EXISTS `org_is_species` (
  `species_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  PRIMARY KEY (`species_id`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `primary_metabolism`
--

DROP TABLE IF EXISTS `primary_metabolism`;
CREATE TABLE IF NOT EXISTS `primary_metabolism` (
  `react_name` varchar(10) NOT NULL,
  `react_id` varchar(10) NOT NULL,
  `react_path` varchar(100) NOT NULL,
  `reacti_ec` varchar(20) NOT NULL,
  `react_ecname` varchar(100) NOT NULL,
  `react_gene` varchar(20) NOT NULL,
  `react_ref` text NOT NULL,
  `react_cbs51388` text NOT NULL,
  `react_atcc1015` text NOT NULL,
  `react_isoenzymes` varchar(20) NOT NULL,
  `react_modelkey` text NOT NULL,
  PRIMARY KEY (`react_name`,`react_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `primary_metabolism_models`
--

DROP TABLE IF EXISTS `primary_metabolism_models`;
CREATE TABLE IF NOT EXISTS `primary_metabolism_models` (
  `model_abb` varchar(100) NOT NULL,
  `model_name` text NOT NULL,
  `model_keggsyn` text NOT NULL,
  `model_keggnr` varchar(10) NOT NULL,
  `model_origin` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `proteins`
--

DROP TABLE IF EXISTS `proteins`;
CREATE TABLE IF NOT EXISTS `proteins` (
  `org_id` int(11) NOT NULL,
  `prot_orgkey` varchar(100) NOT NULL,
  `prot_seqkey` varchar(100) NOT NULL,
  `prot_tailkey` varchar(100) NOT NULL,
  `prot_seqname` varchar(100) NOT NULL,
  `prot_orgseq` varchar(100) DEFAULT NULL,
  `prot_seq` longblob,
  PRIMARY KEY (`prot_seqname`,`org_id`),
  KEY `org_id` (`org_id`),
  KEY `iprot_seq_name` (`prot_seqname`),
  KEY `iprot_seqkey` (`prot_seqkey`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_go`
--

DROP TABLE IF EXISTS `protein_has_go`;
CREATE TABLE IF NOT EXISTS `protein_has_go` (
  `protein_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  `go_term_id` int(11) NOT NULL,
  PRIMARY KEY (`org_id`,`protein_id`,`go_term_id`),
  KEY `goterm_id` (`go_term_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_ipr`
--

DROP TABLE IF EXISTS `protein_has_ipr`;
CREATE TABLE IF NOT EXISTS `protein_has_ipr` (
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `ipr_id` varchar(50) NOT NULL,
  `ipr_domain_start` int(11) NOT NULL,
  `ipr_domain_end` int(11) NOT NULL,
  `ipr_score` varchar(20) NOT NULL,
  PRIMARY KEY (`org_id`,`protein_id`,`ipr_domain_start`,`ipr_domain_end`,`ipr_score`),
  KEY `org_id` (`org_id`),
  KEY `ipr_id` (`ipr_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_kegg`
--

DROP TABLE IF EXISTS `protein_has_kegg`;
CREATE TABLE IF NOT EXISTS `protein_has_kegg` (
  `kegg_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `kegg_ecNum` varchar(50) NOT NULL,
  `kegg_pathway` varchar(500) NOT NULL,
  PRIMARY KEY (`kegg_id`,`protein_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_kog`
--

DROP TABLE IF EXISTS `protein_has_kog`;
CREATE TABLE IF NOT EXISTS `protein_has_kog` (
  `kog_id` varchar(50) NOT NULL,
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  PRIMARY KEY (`kog_id`,`protein_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `sigp`
--

DROP TABLE IF EXISTS `sigp`;
CREATE TABLE IF NOT EXISTS `sigp` (
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `sigp_nn_cutpos` int(11) NOT NULL DEFAULT '0',
  `sigp_neuro_net_vote` int(11) DEFAULT NULL,
  `sigp_hmm_cutpos` int(11) DEFAULT NULL,
  `sigp_hmm_signalpep_probability` float DEFAULT NULL,
  PRIMARY KEY (`org_id`,`protein_id`,`sigp_nn_cutpos`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
CREATE TABLE IF NOT EXISTS `species` (
  `species_id` int(11) NOT NULL AUTO_INCREMENT,
  `species_section` int(11) NOT NULL,
  `species_clade` int(11) NOT NULL,
  `species_desc` varchar(100) NOT NULL,
  PRIMARY KEY (`species_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `testgff`
--

DROP TABLE IF EXISTS `testgff`;
CREATE TABLE IF NOT EXISTS `testgff` (
  `gff_seqorigin` varchar(100) NOT NULL,
  `org_id` int(11) NOT NULL DEFAULT '0',
  `gff_name` varchar(100) NOT NULL,
  `gff_type` varchar(100) NOT NULL,
  `gff_protein_id` varchar(100) DEFAULT NULL,
  `gff_trans_id` varchar(100) DEFAULT NULL,
  `gff_attributes` varchar(500) NOT NULL DEFAULT '',
  `gff_start` int(11) NOT NULL DEFAULT '0',
  `gff_end` int(11) DEFAULT NULL,
  `gff_length` int(11) DEFAULT NULL,
  `gff_score` varchar(10) DEFAULT NULL,
  `gff_strand` varchar(10) DEFAULT NULL,
  `gff_phase` varchar(10) DEFAULT NULL,
  `gff_exonnr` varchar(500) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `tmp`
--

DROP TABLE IF EXISTS `tmp`;
CREATE TABLE IF NOT EXISTS `tmp` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `tmp_rep`
--

DROP TABLE IF EXISTS `tmp_rep`;
CREATE TABLE IF NOT EXISTS `tmp_rep` (
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `qcov` decimal(18,4) DEFAULT NULL,
  `scov` decimal(18,4) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `transcripts`
--

DROP TABLE IF EXISTS `transcripts`;
CREATE TABLE IF NOT EXISTS `transcripts` (
  `org_id` int(11) NOT NULL DEFAULT '0',
  `trans_orgseq` longblob,
  `trans_orgkey` varchar(100) NOT NULL,
  `trans_seqkey` varchar(100) NOT NULL,
  `trans_tailkey` varchar(100) NOT NULL,
  `trans_seqname` varchar(100) NOT NULL,
  `trans_seq` longblob,
  PRIMARY KEY (`trans_seqname`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `assembly`
--
ALTER TABLE `assembly`
  ADD CONSTRAINT `assembly_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `cds`
--
ALTER TABLE `cds`
  ADD CONSTRAINT `cds_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `gff`
--
ALTER TABLE `gff`
  ADD CONSTRAINT `gff_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `org_has_environment`
--
ALTER TABLE `org_has_environment`
  ADD CONSTRAINT `org_has_environment_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON UPDATE CASCADE,
  ADD CONSTRAINT `org_has_environment_ibfk_2` FOREIGN KEY (`env_id`) REFERENCES `environment` (`env_id`) ON UPDATE CASCADE;

--
-- Constraints for table `org_is_species`
--
ALTER TABLE `org_is_species`
  ADD CONSTRAINT `org_is_species_ibfk_3` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `org_is_species_ibfk_5` FOREIGN KEY (`species_id`) REFERENCES `species` (`species_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `proteins`
--
ALTER TABLE `proteins`
  ADD CONSTRAINT `proteins_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `protein_has_go`
--
ALTER TABLE `protein_has_go`
  ADD CONSTRAINT `protein_has_go_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `protein_has_ipr`
--
ALTER TABLE `protein_has_ipr`
  ADD CONSTRAINT `protein_has_ipr_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE,
  ADD CONSTRAINT `protein_has_ipr_ibfk_2` FOREIGN KEY (`ipr_id`) REFERENCES `ipr` (`ipr_id`) ON DELETE CASCADE;

--
-- Constraints for table `protein_has_kegg`
--
ALTER TABLE `protein_has_kegg`
  ADD CONSTRAINT `protein_has_kegg_ibfk_4` FOREIGN KEY (`kegg_id`) REFERENCES `kegg` (`kegg_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `protein_has_kegg_ibfk_5` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `protein_has_kog`
--
ALTER TABLE `protein_has_kog`
  ADD CONSTRAINT `protein_has_kog_ibfk_5` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `protein_has_kog_ibfk_6` FOREIGN KEY (`kog_id`) REFERENCES `kog` (`kog_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `sigp`
--
ALTER TABLE `sigp`
  ADD CONSTRAINT `sigp_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `transcripts`
--
ALTER TABLE `transcripts`
  ADD CONSTRAINT `transcripts_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
