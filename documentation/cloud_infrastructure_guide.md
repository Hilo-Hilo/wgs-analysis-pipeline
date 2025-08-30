# Google Cloud Infrastructure Guide for Genomics Analysis
**Complete Setup and Management Documentation**  
*Updated: July 22, 2025*

## Overview

This guide provides comprehensive instructions for setting up and managing Google Cloud infrastructure for whole genome sequencing analysis, including VM configuration, storage optimization, and data transfer strategies.

## Project Requirements

- **Computational Needs**: 32+ cores, 128GB+ RAM for variant annotation
- **Storage Requirements**: 1TB+ for databases and intermediate files
- **Network**: High bandwidth for database downloads (620GB gnomAD)
- **Budget**: ~$60 total project cost with optimization strategies

## Google Cloud Platform Setup

### Initial Project Configuration

#### 1. Project Creation and Billing
```bash
# Create new GCP project
gcloud projects create wgs-analysis-project --name="WGS CHM13 Analysis"
gcloud config set project wgs-analysis-project

# Enable required APIs
gcloud services enable compute.googleapis.com
gcloud services enable storage.googleapis.com
gcloud services enable logging.googleapis.com

# Set up billing (required for compute resources)
gcloud billing accounts list
gcloud billing projects link wgs-analysis-project --billing-account=XXXXXX-XXXXXX-XXXXXX
```

#### 2. Authentication and Security
```bash
# Install and initialize gcloud CLI
curl https://sdk.cloud.google.com | bash
exec -l $SHELL
gcloud init

# Create service account for automation
gcloud iam service-accounts create wgs-analysis-sa \
  --display-name="WGS Analysis Service Account"

# Grant necessary permissions
gcloud projects add-iam-policy-binding wgs-analysis-project \
  --member="serviceAccount:wgs-analysis-sa@wgs-analysis-project.iam.gserviceaccount.com" \
  --role="roles/compute.admin"
```

### Compute Engine VM Configuration

#### High-Performance VM Specification

**Selected Configuration:**
- **Instance Type**: n2-standard-32
- **vCPUs**: 32 (Intel Cascade Lake)
- **Memory**: 128GB
- **Architecture**: x86_64
- **Network**: Up to 32 Gbps

#### VM Creation Command
```bash
gcloud compute instances create wgs-analysis-vm \
  --zone=us-central1-a \
  --machine-type=n2-standard-32 \
  --network-tier=PREMIUM \
  --maintenance-policy=MIGRATE \
  --provisioning-model=STANDARD \
  --service-account=wgs-analysis-sa@wgs-analysis-project.iam.gserviceaccount.com \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --create-disk=auto-delete=yes,boot=yes,device-name=wgs-analysis-vm,image=projects/ubuntu-os-cloud/global/images/ubuntu-2004-focal-v20240731,mode=rw,size=50,type=projects/wgs-analysis-project/zones/us-central1-a/diskTypes/pd-standard \
  --create-disk=device-name=genomics-data,mode=rw,size=1000,type=projects/wgs-analysis-project/zones/us-central1-a/diskTypes/pd-standard,auto-delete=yes \
  --no-shielded-secure-boot \
  --shielded-vtpm \
  --shielded-integrity-monitoring \
  --labels=project=wgs-analysis,purpose=genomics-compute \
  --reservation-affinity=any
```

**Cost Analysis:**
- **Compute**: $1.76/hour (32 vCPUs × $0.055/hour)
- **Storage**: $40/month (1TB × $0.04/GB/month)
- **Network**: $0.12/GB egress (minimal for this use case)
- **Total Runtime Cost**: ~$14 for 8-hour analysis

### Storage Configuration

#### Persistent Disk Setup

**Disk Specifications:**
- **Primary Boot Disk**: 50GB SSD (Ubuntu 20.04 LTS)
- **Data Disk**: 1TB Standard Persistent Disk
- **Performance**: 3,000 IOPS, 48MB/s throughput

#### Disk Mounting and Formatting
```bash
# SSH into the VM
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a

# Format and mount the data disk
sudo mkfs.ext4 -F /dev/sdb
sudo mkdir -p /mnt/genomics
sudo mount /dev/sdb /mnt/genomics
sudo chown -R $USER:$USER /mnt/genomics

# Make mount persistent
echo '/dev/sdb /mnt/genomics ext4 defaults 0 2' | sudo tee -a /etc/fstab

# Verify mount
df -h /mnt/genomics
# Expected output: /dev/sdb  985G   77M  935G   1% /mnt/genomics
```

#### Storage Expansion Strategy

**Dynamic Disk Resizing:**
```bash
# Expand disk size without downtime
gcloud compute disks resize genomics-data --size=1500GB --zone=us-central1-a

# Resize filesystem on the VM
sudo resize2fs /dev/sdb

# Verify new size
df -h /mnt/genomics
```

**Storage Monitoring:**
```bash
# Monitor disk usage during pipeline
watch -n 30 'df -h /mnt/genomics && du -sh /mnt/genomics/*/'

# Set up alerts for >90% usage
gcloud alpha monitoring policies create \
  --policy-from-file=disk-usage-alert.yaml
```

### Network Configuration and Optimization

#### Firewall Rules
```bash
# Create firewall rule for SSH access
gcloud compute firewall-rules create allow-ssh-wgs \
  --allow tcp:22 \
  --source-ranges 0.0.0.0/0 \
  --description "Allow SSH access to WGS analysis VM" \
  --target-tags wgs-analysis

# Apply firewall tag to VM
gcloud compute instances add-tags wgs-analysis-vm \
  --tags wgs-analysis \
  --zone us-central1-a
```

#### Network Performance Optimization
```bash
# Enable high-bandwidth networking
sudo apt update
sudo apt install -y google-cloud-ops-agent

# Configure network optimizations
echo 'net.core.rmem_default = 262144' | sudo tee -a /etc/sysctl.conf
echo 'net.core.rmem_max = 16777216' | sudo tee -a /etc/sysctl.conf
echo 'net.core.wmem_default = 262144' | sudo tee -a /etc/sysctl.conf
echo 'net.core.wmem_max = 16777216' | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
```

## Data Transfer and Management

### Cloud Storage Integration

#### Storage Bucket Setup
```bash
# Create regional bucket for optimal performance
gsutil mb -p wgs-analysis-project -c STANDARD -l us-central1 gs://wgs-analysis-data

# Set lifecycle policy for cost optimization
cat > lifecycle.json << EOF
{
  "rule": [
    {
      "action": {"type": "SetStorageClass", "storageClass": "NEARLINE"},
      "condition": {"age": 30}
    },
    {
      "action": {"type": "SetStorageClass", "storageClass": "COLDLINE"},
      "condition": {"age": 90}
    }
  ]
}
EOF

gsutil lifecycle set lifecycle.json gs://wgs-analysis-data
```

#### Efficient Data Transfer

**Large File Transfer (Parallel):**
```bash
# Enable parallel composite uploads
gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp \
  bcftools_variants_high_quality.vcf.gz \
  gs://wgs-analysis-data/input/

# Verify transfer integrity
gsutil hash gs://wgs-analysis-data/input/bcftools_variants_high_quality.vcf.gz
md5sum bcftools_variants_high_quality.vcf.gz  # Compare locally
```

**Database Download Optimization:**
```bash
# Use aria2c for parallel downloads
sudo apt install aria2

# Download gnomAD in parallel chunks
aria2c -x 8 -s 8 \
  https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz

# Monitor download progress
watch -n 5 'ls -lah /mnt/genomics/downloads/ && df -h /mnt/genomics'
```

### Data Backup and Versioning

#### Critical Data Backup Strategy
```bash
# Backup essential results
gsutil -m cp /mnt/genomics/results/*.vcf.gz gs://wgs-analysis-data/backup/$(date +%Y%m%d)/
gsutil -m cp /mnt/genomics/results/*.html gs://wgs-analysis-data/backup/$(date +%Y%m%d)/

# Create versioned snapshots of VM disk
gcloud compute snapshots create wgs-analysis-snapshot-$(date +%Y%m%d) \
  --source-disk genomics-data \
  --zone us-central1-a \
  --description "WGS analysis checkpoint $(date +%Y-%m-%d)"
```

#### Automated Backup Script
```bash
#!/bin/bash
# backup_genomics_data.sh - Automated backup script

DATE=$(date +%Y%m%d_%H%M%S)
BACKUP_BUCKET="gs://wgs-analysis-data/backup"

# Backup results
echo "Backing up analysis results..."
gsutil -m cp -r /mnt/genomics/results "${BACKUP_BUCKET}/${DATE}/"

# Create disk snapshot
echo "Creating disk snapshot..."
gcloud compute snapshots create "wgs-analysis-${DATE}" \
  --source-disk genomics-data \
  --zone us-central1-a

echo "Backup completed: ${DATE}"
```

## Resource Monitoring and Management

### System Performance Monitoring

#### Real-Time Resource Tracking
```bash
# Install monitoring tools
sudo apt install -y htop iotop nethogs ncdu

# Create monitoring dashboard script
cat > monitor_resources.sh << 'EOF'
#!/bin/bash
while true; do
  clear
  echo "=== WGS Analysis VM Resource Monitor ==="
  echo "Time: $(date)"
  echo
  echo "=== CPU and Memory ==="
  top -bn1 | head -5
  echo
  echo "=== Disk Usage ==="
  df -h /mnt/genomics
  echo
  echo "=== Active Processes ==="
  ps aux | grep -E 'vep|bcftools|bwa' | grep -v grep
  echo
  echo "=== Network ==="
  ss -tuln | grep :22  # SSH connections
  sleep 30
done
EOF

chmod +x monitor_resources.sh
```

#### Cloud Monitoring Integration
```bash
# Install Cloud Ops Agent
curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
sudo bash add-google-cloud-ops-agent-repo.sh --also-install

# Configure custom metrics
sudo tee /etc/google-cloud-ops-agent/config.yaml > /dev/null << 'EOF'
logging:
  receivers:
    wgs_logs:
      type: files
      include_paths:
        - /mnt/genomics/logs/*.log
      exclude_paths:
        - /mnt/genomics/logs/*.tmp
  service:
    pipelines:
      default_pipeline:
        receivers: [wgs_logs]
        processors: []
        exporters: [google_cloud_logging]

metrics:
  receivers:
    genomics_disk:
      type: disk
      include_devices:
        - /dev/sdb
  service:
    pipelines:
      genomics_pipeline:
        receivers: [genomics_disk]
        processors: []
        exporters: [google_cloud_monitoring]
EOF

sudo systemctl restart google-cloud-ops-agent
```

### Cost Optimization Strategies

#### Preemptible Instance Strategy (Optional)
```bash
# Create cost-optimized preemptible instance (60-91% savings)
gcloud compute instances create wgs-analysis-vm-preempt \
  --zone=us-central1-a \
  --machine-type=n2-standard-32 \
  --preemptible \
  --maintenance-policy=TERMINATE \
  --boot-disk-size=50GB \
  --boot-disk-type=pd-standard

# Note: Use only for fault-tolerant workloads
# Regular instance recommended for critical analysis
```

#### Scheduled Instance Management
```bash
# Stop VM automatically after analysis
# Add to crontab for scheduled shutdown
echo "0 2 * * * gcloud compute instances stop wgs-analysis-vm --zone=us-central1-a" | crontab -

# Start VM on demand
gcloud compute instances start wgs-analysis-vm --zone=us-central1-a
```

#### Storage Cost Optimization
```bash
# Move completed results to cheaper storage classes
gsutil rewrite -s NEARLINE gs://wgs-analysis-data/backup/**
gsutil rewrite -s COLDLINE gs://wgs-analysis-data/archive/**

# Delete temporary files after analysis
find /mnt/genomics/temp -type f -mtime +7 -delete
find /mnt/genomics/downloads -name "*.tmp" -delete
```

## Security and Access Management

### IAM Configuration

#### Principle of Least Privilege
```bash
# Create custom role for genomics analysis
gcloud iam roles create genomicsAnalyst \
  --project=wgs-analysis-project \
  --title="Genomics Analyst" \
  --description="Role for WGS analysis operations" \
  --permissions=compute.instances.get,compute.instances.start,compute.instances.stop,storage.objects.create,storage.objects.delete,storage.objects.get

# Assign role to service account
gcloud projects add-iam-policy-binding wgs-analysis-project \
  --member="serviceAccount:wgs-analysis-sa@wgs-analysis-project.iam.gserviceaccount.com" \
  --role="projects/wgs-analysis-project/roles/genomicsAnalyst"
```

#### SSH Key Management
```bash
# Generate project-specific SSH keys
ssh-keygen -t rsa -b 4096 -f ~/.ssh/wgs-analysis-key -C "wgs-analysis@project"

# Add public key to project metadata
gcloud compute project-info add-metadata \
  --metadata-from-file ssh-keys=<(echo "$USER:$(cat ~/.ssh/wgs-analysis-key.pub)")

# Connect using specific key
ssh -i ~/.ssh/wgs-analysis-key $USER@$(gcloud compute instances describe wgs-analysis-vm --zone=us-central1-a --format="value(networkInterfaces[0].accessConfigs[0].natIP)")
```

### Data Security

#### Disk Encryption
```bash
# Create encrypted disk
gcloud compute disks create genomics-data-encrypted \
  --size=1000GB \
  --zone=us-central1-a \
  --type=pd-standard \
  --kms-key=projects/wgs-analysis-project/locations/us-central1/keyRings/genomics-ring/cryptoKeys/genomics-key

# Attach encrypted disk to VM
gcloud compute instances attach-disk wgs-analysis-vm \
  --disk=genomics-data-encrypted \
  --zone=us-central1-a
```

#### Network Security
```bash
# Create VPC with restricted access
gcloud compute networks create wgs-vpc --subnet-mode=custom

gcloud compute networks subnets create wgs-subnet \
  --network=wgs-vpc \
  --range=10.0.0.0/24 \
  --region=us-central1

# Firewall rules for restricted access
gcloud compute firewall-rules create wgs-allow-internal \
  --network=wgs-vpc \
  --allow=tcp,udp,icmp \
  --source-ranges=10.0.0.0/24

gcloud compute firewall-rules create wgs-allow-ssh \
  --network=wgs-vpc \
  --allow=tcp:22 \
  --source-ranges=YOUR_IP_ADDRESS/32  # Replace with your IP
```

## Troubleshooting and Common Issues

### VM Connection Issues

#### SSH Connection Problems
```bash
# Debug SSH connectivity
gcloud compute instances list
gcloud compute firewall-rules list --filter="name:ssh"

# Reset SSH keys if needed
gcloud compute reset-windows-password wgs-analysis-vm --zone=us-central1-a --user=$USER

# Use browser-based SSH if local connection fails
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --ssh-flag="-v"
```

#### VM Performance Issues
```bash
# Check VM status and health
gcloud compute instances describe wgs-analysis-vm --zone=us-central1-a

# Monitor resource usage
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="top -bn1"

# Check for system logs
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="sudo journalctl -u google-startup-scripts.service"
```

### Storage Issues

#### Disk Space Problems
```bash
# Check disk usage
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="df -h"

# Identify large files
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="du -h /mnt/genomics | sort -rh | head -20"

# Clean up temporary files
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="find /tmp -type f -mtime +1 -delete"
```

#### Disk Mounting Issues
```bash
# Check disk attachment
gcloud compute instances describe wgs-analysis-vm --zone=us-central1-a --format="table(disks[].deviceName,disks[].source)"

# Remount disk if needed
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="sudo umount /mnt/genomics && sudo mount /dev/sdb /mnt/genomics"

# Verify filesystem integrity
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="sudo fsck -f /dev/sdb"
```

### Network and Download Issues

#### Slow Download Speeds
```bash
# Test network performance
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="wget -O /dev/null http://speedtest.tele2.net/100MB.zip"

# Use multiple connections for downloads
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="aria2c -x 8 -s 8 [URL]"

# Check DNS resolution
gcloud compute ssh wgs-analysis-vm --zone=us-central1-a --command="nslookup gnomad-public-us-east-1.s3.amazonaws.com"
```

## Project Cleanup and Cost Management

### Resource Cleanup

#### Complete Project Shutdown
```bash
# Stop all running instances
gcloud compute instances stop --all --zones=us-central1-a

# Delete instances (WARNING: This removes all data)
gcloud compute instances delete wgs-analysis-vm --zone=us-central1-a

# Delete disks
gcloud compute disks delete genomics-data --zone=us-central1-a

# Delete snapshots older than 30 days
gcloud compute snapshots list --filter="creationTimestamp<$(date -d '30 days ago' -Iseconds)" --format="value(name)" | xargs -r gcloud compute snapshots delete
```

#### Selective Resource Management
```bash
# Stop VM but keep disks
gcloud compute instances stop wgs-analysis-vm --zone=us-central1-a

# Delete temporary storage buckets
gsutil rm -r gs://wgs-analysis-data/temp/**

# Archive completed results
gsutil mv gs://wgs-analysis-data/results/** gs://wgs-analysis-data/archive/
```

### Cost Monitoring

#### Billing Alerts
```bash
# Set up budget alert
gcloud billing budgets create \
  --billing-account=XXXXXX-XXXXXX-XXXXXX \
  --display-name="WGS Analysis Budget" \
  --budget-amount=100USD \
  --threshold-rule=percent=50,basis=CURRENT_SPEND \
  --threshold-rule=percent=90,basis=CURRENT_SPEND
```

#### Cost Analysis
```bash
# Export billing data
gcloud billing accounts list
gcloud billing export describe --billing-account=XXXXXX-XXXXXX-XXXXXX

# Query costs using BigQuery
bq query --use_legacy_sql=false '
SELECT
  service.description as service,
  SUM(cost) as total_cost
FROM `PROJECT_ID.billing.gcp_billing_export_v1_BILLING_ACCOUNT_ID`
WHERE project.id = "wgs-analysis-project"
GROUP BY service
ORDER BY total_cost DESC;'
```

## Best Practices Summary

### Infrastructure Management
1. **Right-size resources** based on workload requirements
2. **Use persistent disks** for data that must survive VM termination
3. **Implement automated backups** for critical results
4. **Monitor costs** with budgets and alerts
5. **Apply security best practices** with IAM and encryption

### Performance Optimization
1. **Choose appropriate regions** close to data sources
2. **Use parallel processing** for large file transfers
3. **Optimize network settings** for genomics workloads
4. **Monitor resource utilization** and adjust accordingly
5. **Clean up temporary files** regularly to save costs

### Cost Management
1. **Stop VMs when not in use** (70-80% savings)
2. **Use preemptible instances** for fault-tolerant workloads
3. **Implement lifecycle policies** for storage optimization
4. **Delete unnecessary snapshots** and backups regularly
5. **Monitor usage patterns** and optimize resource allocation

---

**Document Status**: Production Ready  
**Last Updated**: July 22, 2025  
**Infrastructure Version**: Google Cloud Platform Optimized  
**Next Review**: Monthly cost and performance analysis
