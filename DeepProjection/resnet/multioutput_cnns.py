import torch
import torch.nn as nn
import torchvision.models as models
import torch.nn.functional as F
from resnet import resnet18, resnet34, resnet50

# 3-layer Multi-Output CNN (Late)
class MultiOutputCNN_3Layer(nn.Module):
    def __init__(self, num_particles=11, num_counts=4, hidden_dim=8):
        super(MultiOutputCNN_3Layer, self).__init__()
        self.conv1 = nn.Conv2d(1, hidden_dim, 2, 2) # (8, 64, 64)
        self.conv2 = nn.Conv2d(hidden_dim, hidden_dim * 4, 4, 4) # (32, 16, 16)
        self.conv3 = nn.Conv2d(hidden_dim * 4, hidden_dim * 16, 4, 4) # (128, 4, 4)
        self.dropout1 = nn.Dropout(0.25)
        self.fc1 = nn.Linear(2048, num_counts)
        self.fc2 = nn.Linear(2048, num_particles)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)

        y1 = self.fc1(x)
        y1 = F.log_softmax(y1, dim=1)
        y2 = self.fc2(x)
        y2 = F.log_softmax(y2, dim=1)
        return y1, y2
    
# 5-Layer Multi-Output CNN (Late)
class MultiOutputCNN_5Layer(nn.Module):
    def __init__(self, num_particles=11, num_counts=4, hidden_dim=8):
        super(MultiOutputCNN_5Layer, self).__init__()
        self.conv1 = nn.Conv2d(1, hidden_dim, 2, 2) # (8, 64, 64)
        self.conv2 = nn.Conv2d(hidden_dim, hidden_dim * 2, 2, 2) # (16, 32, 32)
        self.conv3 = nn.Conv2d(hidden_dim * 2, hidden_dim * 4, 2, 2) # (32, 16, 16)
        self.conv4 = nn.Conv2d(hidden_dim * 4, hidden_dim * 8, 2, 2) # (64, 8, 8)
        self.conv5 = nn.Conv2d(hidden_dim * 8, hidden_dim * 16, 2, 2) # (128, 4, 4)
        self.dropout1 = nn.Dropout(0.25)
        self.dropout2 = nn.Dropout(0.5)
        self.fc1 = nn.Linear(2048, num_counts)
        self.fc2 = nn.Linear(2048, num_particles)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)
        x = self.conv4(x)
        x = F.relu(x)
        x = self.conv5(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)

        y1 = self.fc1(x)
        y1 = F.log_softmax(y1, dim=1)
        y2 = self.fc2(x)
        y2 = F.log_softmax(y2, dim=1)
        return y1, y2
    
# 10-Layer Multi-Output CNN (Late)
class MultiOutputCNN_10Layer(nn.Module):
    def __init__(self, num_particles=11, num_counts=4, hidden_dim=8):
        super(MultiOutputCNN_10Layer, self).__init__()
        self.conv1 = nn.Conv2d(1, hidden_dim, 2, 2) # (8, 64, 64)
        self.conv2 = nn.Conv2d(hidden_dim, hidden_dim * 2, 2, 2) # (16, 32, 32)
        self.conv3 = nn.Conv2d(hidden_dim * 2, hidden_dim * 4, 2, 2) # (32, 16, 16)
        self.conv4 = nn.Conv2d(hidden_dim * 4, hidden_dim * 8, 2, 2) # (64, 8, 8)
        self.conv5 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 7, 7)
        self.conv6 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 6, 6)
        self.conv7 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 5, 5)
        self.conv8 = nn.Conv2d(hidden_dim * 8, hidden_dim * 16, 1, 1) # (128, 5, 5)
        self.conv9 = nn.Conv2d(hidden_dim * 16, hidden_dim * 16, 2, 1) # (128, 4, 4)
        self.dropout1 = nn.Dropout(0.25)
        self.fc1 = nn.Linear(2048, num_counts)
        self.fc2 = nn.Linear(2048, num_particles)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)
        x = self.conv4(x)
        x = F.relu(x)
        x = self.conv5(x)
        x = F.relu(x)
        x = self.conv6(x)
        x = F.relu(x)
        x = self.conv7(x)
        x = F.relu(x)
        x = self.conv8(x)
        x = F.relu(x)
        x = self.conv9(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)

        y1 = self.fc1(x)
        y1 = F.log_softmax(y1, dim=1)
        y2 = self.fc2(x)
        y2 = F.log_softmax(y2, dim=1)
        return y1, y2
    
# 18-Layer Multi-Output CNN
class MultiOutputCNN_18Layer(nn.Module):
    def __init__(self, num_particles=11, num_counts=4, hidden_dim=8):
        super(MultiOutputCNN_18Layer, self).__init__()
        self.conv1 = nn.Conv2d(1, hidden_dim, 2, 2) # (8, 64, 64)
        self.conv2 = nn.Conv2d(hidden_dim, hidden_dim * 2, 2, 2) # (16, 32, 32)
        self.conv3 = nn.Conv2d(hidden_dim * 2, hidden_dim * 4, 2, 2) # (32, 16, 16)
        self.conv4 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 15, 15)
        self.conv5 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 14, 14)
        self.conv6 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 13, 13)
        self.conv7 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 12, 12)
        self.conv8 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 11, 11)
        self.conv9 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 10, 10)
        self.conv10 = nn.Conv2d(hidden_dim * 4, hidden_dim * 4, 2, 1) # (32, 9, 9)
        self.conv11 = nn.Conv2d(hidden_dim * 4, hidden_dim * 8, 1, 1) # (64, 9, 9)
        self.conv12 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 8, 8)
        self.conv13 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 7, 7)
        self.conv14 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 6, 6)
        self.conv15 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 5, 5)
        self.conv16 = nn.Conv2d(hidden_dim * 8, hidden_dim * 16, 1, 1) # (128, 5, 5)
        self.conv17 = nn.Conv2d(hidden_dim * 16, hidden_dim * 16, 2, 1) # (128, 4, 4)
        self.dropout1 = nn.Dropout(0.25)
        self.fc1 = nn.Linear(2048, num_counts)
        self.fc2 = nn.Linear(2048, num_particles)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)
        x = self.conv4(x)
        x = F.relu(x)
        x = self.conv5(x)
        x = F.relu(x)
        x = self.conv6(x)
        x = F.relu(x)
        x = self.conv7(x)
        x = F.relu(x)
        x = self.conv8(x)
        x = F.relu(x)
        x = self.conv9(x)
        x = F.relu(x)
        x = self.conv10(x)
        x = F.relu(x)
        x = self.conv11(x)
        x = F.relu(x)
        x = self.conv12(x)
        x = F.relu(x)
        x = self.conv13(x)
        x = F.relu(x)
        x = self.conv14(x)
        x = F.relu(x)
        x = self.conv15(x)
        x = F.relu(x)
        x = self.conv16(x)
        x = F.relu(x)
        x = self.conv17(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)

        y1 = self.fc1(x)
        y1 = F.log_softmax(y1, dim=1)
        y2 = self.fc2(x)
        y2 = F.log_softmax(y2, dim=1)
        return y1, y2
    
# Multi-Output CNN (Early)
class MultiOutputCNN_Early(nn.Module):
    def __init__(self, num_particles=11, num_counts=4, hidden_dim=8):
        super(MultiOutputCNN_Early, self).__init__()
        
        
        self.conv1 = nn.Conv2d(1, hidden_dim, 2, 2) # (8, 64, 64)
        self.conv2 = nn.Conv2d(hidden_dim, hidden_dim * 2, 2, 2) # (16, 32, 32)
        
        #Branched
        self.conv3_b1 = nn.Conv2d(hidden_dim * 2, hidden_dim * 4, 2, 2) # (32, 16, 16)
        self.conv4_b1 = nn.Conv2d(hidden_dim * 4, hidden_dim * 8, 2, 2) # (64, 8, 8)
        self.conv5_b1 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 7, 7)
        self.conv6_b1 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 6, 6)
        self.conv7_b1 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 5, 5)
        self.conv8_b1 = nn.Conv2d(hidden_dim * 8, hidden_dim * 16, 1, 1) # (128, 5, 5)
        self.conv9_b1 = nn.Conv2d(hidden_dim * 16, hidden_dim * 16, 2, 1) # (128, 4, 4)
        self.dropout1_b1 = nn.Dropout(0.25)
        self.fc1_b1 = nn.Linear(2048, num_counts)
        
        self.conv3_b2 = nn.Conv2d(hidden_dim * 2, hidden_dim * 4, 2, 2) # (32, 16, 16)
        self.conv4_b2 = nn.Conv2d(hidden_dim * 4, hidden_dim * 8, 2, 2) # (64, 8, 8)
        self.conv5_b2 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 7, 7)
        self.conv6_b2 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 6, 6)
        self.conv7_b2 = nn.Conv2d(hidden_dim * 8, hidden_dim * 8, 2, 1) # (64, 5, 5)
        self.conv8_b2 = nn.Conv2d(hidden_dim * 8, hidden_dim * 16, 1, 1) # (128, 5, 5)
        self.conv9_b2 = nn.Conv2d(hidden_dim * 16, hidden_dim * 16, 2, 1) # (128, 4, 4)
        self.dropout1_b2 = nn.Dropout(0.25)
        self.fc1_b2 = nn.Linear(2048, num_particles)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        
        #Branched
        y1 = self.conv3_b1(x)
        y1 = F.relu(y1)
        y1 = self.conv4_b1(y1)
        y1 = F.relu(y1)
        y1 = self.conv5_b1(y1)
        y1 = F.relu(y1)
        y1 = self.conv6_b1(y1)
        y1 = F.relu(y1)
        y1 = self.conv7_b1(y1)
        y1 = F.relu(y1)
        y1 = self.conv8_b1(y1)
        y1 = F.relu(y1)
        y1 = self.conv9_b1(y1)
        y1 = F.relu(y1)
        y1 = self.dropout1_b1(y1)
        y1 = torch.flatten(y1, 1)
        y1 = self.fc1_b1(y1)
        y1 = F.log_softmax(y1, dim=1)

        y2 = self.conv3_b2(x)
        y2 = F.relu(y2)
        y2 = self.conv4_b2(y2)
        y2 = F.relu(y2)
        y2 = self.conv5_b2(y2)
        y2 = F.relu(y2)
        y2 = self.conv6_b2(y2)
        y2 = F.relu(y2)
        y2 = self.conv7_b2(y2)
        y2 = F.relu(y2)
        y2 = self.conv8_b2(y2)
        y2 = F.relu(y2)
        y2 = self.conv9_b2(y2)
        y2 = F.relu(y2)
        y2 = self.dropout1_b2(y2)
        y2 = torch.flatten(y2, 1)
        y2 = self.fc1_b2(y2)
        y2 = F.log_softmax(y2, dim=1)

        return y1, y2
    
# Multi-Output ResNet18
class CustomResNet18Model(nn.Module):
    def __init__(self, num_counts, num_particles):
        super(CustomResNet18Model, self).__init__()
        self.model_resnet = models.resnet18(pretrained=False)
        self.model_resnet.conv1 = torch.nn.Conv1d(1, 64, (7, 7), (2, 2), (3, 3), bias=True)
        
        self.model_resnet.fc.register_forward_hook(lambda m, inp, out: F.dropout(out, p=0.5, training=m.training))
        
        num_ftrs = self.model_resnet.fc.in_features
        self.model_resnet.fc = nn.Identity()
        self.fc1 = nn.Linear(num_ftrs, num_counts)
        self.fc2 = nn.Linear(num_ftrs, num_particles)
    def forward(self, x):
        x = self.model_resnet(x)
        out1 = self.fc1(x)
        y1 = F.log_softmax(out1, dim=1)
        out2 = self.fc2(x)
        y2 = F.log_softmax(out2, dim=1)
        return y1, y2

# Multi-Output VGG16
class CustomVgg16Model(nn.Module):
    def __init__(self, num_counts, num_particles):
        super(CustomVgg16Model, self).__init__()
        self.model_vgg16 = models.vgg16(pretrained=False, progress=True)
        self.model_vgg16.features[0] = torch.nn.Conv2d(1, 64, (3, 3), (1, 1), (1, 1))
        num_ftrs = self.model_vgg16.classifier[0].in_features
        self.model_vgg16.classifier = nn.Identity()
        self.fc1 = nn.Linear(num_ftrs, num_counts)
        self.fc2 = nn.Linear(num_ftrs, num_particles)
    def forward(self, x):
        x = self.model_vgg16(x)
        out1 = self.fc1(x)
        y1 = F.log_softmax(out1, dim=1)
        out2 = self.fc2(x)
        y2 = F.log_softmax(out2, dim=1)
        return y1, y2