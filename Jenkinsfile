pipeline {
    options {
        timestamps()
    }
    parameters {
        string(name: 'BUILD_VERSION', defaultValue: '', description: 'The build version to deploy (optional)')
        string(name: 'IMAGE_VERSION', defaultValue: '', description: 'The Docker image version to use (mandatory)')
        string(name: 'ENVIRONMENT', defaultValue: 'ci', description: 'Role Name (mandatory)')
    }
    agent {
        label 'ci && deploy && sctl'
    }
    environment {
        PROJECT_NAME     = "sctl"
        DOCKER_REPO_NAME = "ncats/sctl-rshiny-complex"
        INIT_TOKEN       = credentials('Vault-Access')
        SPHINX_TOKEN     = credentials('ncatssvcdvops-sphinx')
        ROLE_NAME        = "$ENVIRONMENT-$PROJECT_NAME"
        APP_TYPE         = "keys"
        DOCKER_IMAGE     = "rstudio/rstudio-connect:${params.IMAGE_VERSION}"
    }
    stages {
        stage('Docker/Apps getSecrets By Role') {
            steps {
                cleanWs()
                checkout scm
                script {
                    sh '''
                    git clone https://$SPHINX_TOKEN@github.com/Sphinx-Automation/devops-pipeline-artifacts.git
                    cd devops-pipeline-artifacts/application
                    /bin/bash getDockerHubSecretsByRole.sh
                    /bin/bash getAppSecretsByRole.sh
                    '''
                }
            }
        }
        stage('Build Version') {
            steps {
                script {
                    env.BUILD_VERSION = VersionNumber(
                        versionNumberString: 'v${BUILD_YEAR, XX}.${BUILD_MONTH, XX}.${BUILD_DAY, XX}.${BUILDS_TODAY}',
                        projectStartDate:    '1970-01-01',
                        skipFailedBuilds:    true)
                    currentBuild.displayName = env.BUILD_VERSION
                }
            }
        }
        stage('Build') {
            steps {
                configFileProvider([
                    configFile(fileId: 'prepare.sh', targetLocation: 'prepare.sh')
                ]) {
                    script {
                        sh '''#!/bin/bash
                        source prepare.sh
                        docker login -u "${DOCKERLOGIN}" -p "${DOCKERPASSWORD}"
                        docker pull ${DOCKER_IMAGE}
                        docker tag ${DOCKER_IMAGE} ${DOCKER_REPO_NAME}:${BUILD_VERSION}
                        docker push ${DOCKER_REPO_NAME}:${BUILD_VERSION}
                        '''
                    }
                }
            }
        }
        stage('Deploy') {
            steps {
                configFileProvider([
                    configFile(fileId: 'docker-compose.yml', targetLocation: 'docker-compose.yml'),
                    configFile(fileId: 'rstudio-connect.gcfg', targetLocation: 'rstudio-connect.gcfg')
                ]) {
                    script {
                        sh '''#!/bin/bash
                        cd /home/ops/rstudio-connect
                        docker-compose down -v --rmi all
                        docker-compose -p ${PROJECT_NAME} up -d
                        '''
                    }
                }
            }
        }
    }
    post {
        always {
            echo "Clean up the workspace in deploy node!"
            cleanWs()
        }
    }
}